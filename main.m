clc
close all
clear all
%% Room configuration
room.length = 4;
room.width = 4;
room.height = 2.5;

%% LEDs scenario configuration
%% configuration 1
% led.position_x = [1.25, 3.75];
% led.position_y = [1.25, 3.75]';
% led.position_z = 3;

% led.position_matrix_x = [1.25,3.75;1.25,3.75];
% led.position_matrix_y = [3.75,3.75;1.25,1.25];

led.position_x = [1,3,1,3];
led.position_y = [3,3,1,1];

% led.position_x = [0.5,0.5,1.5,1.5];
% led.position_y = [0.5,1.5,0.5,1.5];
%%
theta = 70;     % semi-angle at half power
%m = log10(2)/log10(cosd(theta));    %Lambertian order of emission
m = 20;
Adet = 1e-4;    %detector physical area of a PD
global P_LED;

%% receiver configuration, point, or plane, or objects
% receiver plan
space = 0.25;
[X, Y] = meshgrid(0:space:room.length, 0:space:room.width);
[r,c] = size(X);
receiver_height = 0.85;

Ts=1; %gain of an optical filter; ignore if no filter is used 
index=1.5; %refractive index of a lens at a PD; ignore if no lens is used
FOV=85*pi/180; %FOV of a receiver
G_Con=(index^2)/sin(FOV); %gain of an optical concentrator


%% Figure room
% axis([0 5 0 3]);

[XX,YY,ZZ] = meshgrid(0:room.length);
V = XX.*exp(-XX.^2-YY.^2-ZZ.^2);
xslice = [];%0.8;   
yslice = [];
zslice = receiver_height;
figure('Color','w')
slice(XX,YY,ZZ,V,xslice,yslice,zslice,'nearest');colormap white;
axis([0 room.length 0 room.width 0 room.height]);
hold on
%figure('Color','w')
led.position_z = [room.height,room.height,room.height,room.height];
plot3(led.position_x,led.position_y,led.position_z,'gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','y');
%axis([0 room.length 0 room.width 0 room.height]);
grid on
LEDs={'LED1','LED2','LED3','LED4'};
text(led.position_x+0.15,led.position_y+0.15,led.position_z,LEDs);
% hold on
% plot3(X(9,9),Y(9,9),receiver_height,'gs',...
%     'LineWidth',2,...
%     'MarkerSize',10,...
%     'MarkerEdgeColor','k',...
%     'MarkerFaceColor','g');
% PD={'receiver'};
% text(X(9,9)+0.3,Y(9,9)+0.3,receiver_height,PD);
hold off

%% ================   signal generation  ==================================
%% ========         Channel Definition   ==================================
% Create a LED model based on a low-pass filter.
% Indicate the signal bandwidth and measured LED cutoff frequency.
signalBandwidth = 20*10^6;
ledCutoffFrequecy = 3.5*10^6;

% x = ledCutoffFrequecy/signalBandwidth;

ledFilterOrder = 31;
t = 0:ledFilterOrder;
led_filter = exp(-2*pi*(ledCutoffFrequecy/(8*10^6))*t');
led_filter = led_filter/sum(led_filter);
% fvtool(led_filter,1);

%% Initialize the transceiver
% Create a QAM modulator System object(TM) and set the modulation order to
% 16.
numberOfSubcarrier = 4;
baseOversampling = 5;
if numberOfSubcarrier == 1
    %     M = 64; % ori
    M = 16;
    M_psk = 2;
    numSamplesPerSymbol = baseOversampling*numberOfSubcarrier; % numSamplesPerSymbol = up-sampling factor ns
elseif numberOfSubcarrier == 2
    M = [16 4];
    M_psk = [2 2];
    numSamplesPerSymbol = baseOversampling*numberOfSubcarrier;
elseif numberOfSubcarrier == 3
    M = [16 16 4];
    M_psk = [2 2 2];
    numSamplesPerSymbol = baseOversampling*numberOfSubcarrier;
elseif numberOfSubcarrier == 4
    %M = [128 64 32 4];
    M = [128 64 32 4];
    M_psk = [2 2 2 2];
    numSamplesPerSymbol = baseOversampling*numberOfSubcarrier;
elseif numberOfSubcarrier == 5
    %M = [128 64 32 4];
    M = [128 64 32 4];
    M_psk = [2 2 2 2 2];
    numSamplesPerSymbol = baseOversampling*numberOfSubcarrier;
elseif numberOfSubcarrier == 8
    M = [32 16 4 16 16 16 4 4];
    numSamplesPerSymbol = baseOversampling*numberOfSubcarrier;
elseif numberOfSubcarrier == 10
    M = [128 64 32 16 4 4 4 4 2 2];
    numSamplesPerSymbol = baseOversampling*numberOfSubcarrier;
elseif numberOfSubcarrier == 20
    M = [256 256 128 128 128 64 64 64 64 32 32 16 16 16 4 4 4 4 4 2];
    % M = [64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64 64];
    numSamplesPerSymbol = baseOversampling*numberOfSubcarrier;
end
%k = log2(M);
k = log2(M_psk);
n = 6000;
E_bN_0=20;   % Define range of Eb/N0 for calculations
SNR = E_bN_0 + 10*log10(log2(2)); %Conversion into SNR
%SNR = 30;
% % % IQ Imbalance compensator
% % hIQComp = comm.IQImbalanceCompensator('StepSizeSource','Input port', ...
% %             'CoefficientOutputPort',true);
% % stepSize = 1e-5;

%% Define the Raised Cosine Filter for CAP modulation
rolloff = 0.5; %0.9;
span = 10; % I/Q filters length
delay = span;
filterTaps = span*numSamplesPerSymbol+1;

%% Generate Data
% Create a column vector of M-ary random integer symbols.
data = zeros(n/numberOfSubcarrier,numberOfSubcarrier);
modData = zeros(n/numberOfSubcarrier,numberOfSubcarrier);
txSignal = zeros((n/numberOfSubcarrier)*numSamplesPerSymbol,numberOfSubcarrier);
txSignalPe = zeros((n/numberOfSubcarrier)*numSamplesPerSymbol,numberOfSubcarrier);
I_filter = zeros(numberOfSubcarrier,filterTaps);
Q_filter = zeros(numberOfSubcarrier,filterTaps);
inphaseGain = zeros(numberOfSubcarrier);
quadGain = zeros(numberOfSubcarrier,1);


for i = 1:numberOfSubcarrier
    
    %% Create data for each subcarrier
    %     rng(i);
    %     if i==1
    %         data_pre(:,i) = randi([0 M(i)-1],(n/numberOfSubcarrier)-10,1);
    %     elseif i==2
    %         %data_pre(:,i) = zeros((n/numberOfSubcarrier)-10,1);
    %         data_pre(:,i) = randi([0 M(i)-1],(n/numberOfSubcarrier)-10,1);
    %     elseif i==3
    %         %data_pre(:,i) = zeros((n/numberOfSubcarrier)-10,1);
    %         data_pre(:,i) = randi([0 M(i)-1],(n/numberOfSubcarrier)-10,1);
    %     elseif i==4
    %         data_pre(:,i) = randi([0 M(i)-1],(n/numberOfSubcarrier)-10,1);
    %     end
    %     SF = ones(10,1);
    %     data(:,i) = vertcat(SF,data_pre(:,i));
    
    % rng(i);
    if i==1
        data_pre(:,i) = randi([0 M_psk(i)-1],(n/numberOfSubcarrier),1);
    elseif i==2
        %data_pre(:,i) = ones((n/numberOfSubcarrier),1);
        data_pre(:,i) = randi([0 M_psk(i)-1],(n/numberOfSubcarrier),1);
    elseif i==3
        %data_pre(:,i) = ones((n/numberOfSubcarrier),1);
        data_pre(:,i) = randi([0 M_psk(i)-1],(n/numberOfSubcarrier),1);
    elseif i==4
        data_pre(:,i) = randi([0 M_psk(i)-1],(n/numberOfSubcarrier),1);
    elseif i==5
        data_pre(:,i) = randi([0 M_psk(i)-1],(n/numberOfSubcarrier),1);
    end
    
    data(:,i) = data_pre(:,i);
    
    
    %% CAP modulation
    
    %modData(:,i) = qammod(data(:,i),M(i));
    modData(:,i) = pskmod(data(:,i),M_psk(i));
    [I_filter(i,:),Q_filter(i,:),inphaseGain(i),quadGain(i)] = capcosinefilter(rolloff,span,numSamplesPerSymbol,i);
    [txSignal(:,i)] = capmod(modData(:,i),numSamplesPerSymbol,I_filter(i,:),Q_filter(i,:));
    
end


if numberOfSubcarrier == 1
    snr_sc(1) = 26;
    for i = 1:numberOfSubcarrier
        txSignal(:,i) = awgn(txSignal(:,i), snr_sc(i), 'measured');
    end
    
elseif numberOfSubcarrier == 10
    snr_sc(1) = 23;
    snr_sc(2) = 17;
    snr_sc(3) = 13;
    snr_sc(4) = 12;
    snr_sc(5) = 8;
    snr_sc(6) = 7;
    snr_sc(7) = 6;
    snr_sc(8) = 3;
    snr_sc(9) = 2;
    snr_sc(10) = 0;
    for i = 1:numberOfSubcarrier
        txSignal(:,i) = awgn(txSignal(:,i), snr_sc(i), 'measured');
    end
elseif numberOfSubcarrier == 20
    snr_sc(1) = 35;%
    snr_sc(2) = 33;
    snr_sc(3) = 32;%
    snr_sc(4) = 31;
    snr_sc(5) = 29;%
    snr_sc(6) = 26;
    snr_sc(7) = 24;%
    snr_sc(8) = 23;
    snr_sc(9) = 21;%
    snr_sc(10) = 19;
    snr_sc(11) = 17;%
    snr_sc(12) = 14;
    snr_sc(13) = 11;%
    snr_sc(14) = 9;
    snr_sc(15) = 7;%
    snr_sc(16) = 7;
    snr_sc(17) = 6;%
    snr_sc(18) = 6;
    snr_sc(19) = 5;%
    snr_sc(20) = 5;
    for i = 1:numberOfSubcarrier
        txSignal(:,i) = awgn(txSignal(:,i), snr_sc(i), 'measured');
    end
    % elseif numberOfSubcarrier == 4
    %     snr_sc(1) = 35;
    %     snr_sc(2) = 26;
    %     snr_sc(3) = 17;
    %     snr_sc(4) = 7;
    %     for i = 1:numberOfSubcarrier
    %         txSignal(:,i) = awgn(txSignal(:,i), snr_sc(i), 'measured');
    %     end
end

%% Transmitted Signal
txSignal = txSignal';
if numberOfSubcarrier > 1
    tx = sum(txSignal);
else
    tx = txSignal';
end

%% Calculate Spectrum of transmitted signal
Fs = signalBandwidth;       % Sampling frequency
T = 1/Fs;               % Sampling period
L = length(tx);         % Length of signal
time = (0:L-1)*T;       % Time vector

nfft = 256;
[pxx, freq] = pwelch(tx(1,:),hanning(nfft),nfft/2,nfft,Fs);
%P_LED = [pxx(7),pxx(20),pxx(33),pxx(46)];
P_LED = [pxx(14),pxx(40),pxx(66),pxx(72)];
figure(3);
set(gca,'Color','w');
plot(freq,smooth(10*log10(pxx)),'linewidth',1.2);
grid on
title('PSD of transmitted signal')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
% xlim([0 10*1E6])
% ylim([-150 -40])

%% LED coverage
radius = tand(theta)*(room.height-receiver_height);

%%%%%%%%%%%%%%%%%%%%% Plot coverage %%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w');
rectangle('Position',[0,0,room.length,room.width],'FaceColor',[0.85 0.85 0.85],'EdgeColor','b',...
    'LineWidth',3);%  [x,y,w,h] 给定起点[x,y]  矩形宽w高h
hold on

plot(led.position_x,led.position_y,'gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','y');
LEDs={'LED1','LED2','LED3','LED4'};
text(led.position_x(1,:)+0.15,led.position_y(1,:),LEDs);
hold on

[X1,Y1]=circle([led.position_x(1),led.position_y(1)],radius,150);
plot(X1,Y1,'+')
xlabel('Length of surface [m]')
ylabel('Breadth of surface [m]')
title('Top view of LED light coverage')
hold on

[X1,Y1]=circle([led.position_x(2),led.position_y(2)],radius,150);
plot(X1,Y1,'+')
hold on

[X1,Y1]=circle([led.position_x(3),led.position_y(3)],radius,150);
plot(X1,Y1,'+')
hold on

[X1,Y1]=circle([led.position_x(4),led.position_y(4)],radius,150);
plot(X1,Y1,'+')
hold off
legend('LEDs')

%%%%%%%%%%%%%%%%%%%%% End Circle %%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculate Pr for each LED
nLED = 4;

Pt = P_LED;
Pr = zeros(r,c,nLED);
Hlos_fromEachLed = zeros(r,c,nLED);
Pr_led1 = zeros(r,c);
Pr_led2 = zeros(r,c);
Pr_led3 = zeros(r,c);
Pr_led4 = zeros(r,c);
Hlos_led1 = zeros(r,c);
Hlos_led2 = zeros(r,c);
Hlos_led3 = zeros(r,c);
Hlos_led4 = zeros(r,c);
locx = zeros(r,c);
locy = zeros(r,c);
locz = zeros(r,c);
error = zeros(r,c);

for i = 1:r
    for j = 1:c
        for index_led = 1:nLED
            d = sqrt((X(i,j)-led.position_x(index_led))^2 + (Y(i,j)-led.position_y(index_led))^2 + (receiver_height-room.height)^2);
            cosphi = (room.height-receiver_height)/d;
            
            Hlos = (m+1)*Adet*cosphi^(m+1)/(2*pi*d^2);
            Hlos_fromEachLed(i,j,index_led) = Hlos;
            %Pr(i,j,index_led) = Pt(index_led)*Pr_LED(index_led)*Hlos;
        end

        Hlos_led1(i,j) = Hlos_fromEachLed(i,j,1);
        Hlos_led2(i,j) = Hlos_fromEachLed(i,j,2);
        Hlos_led3(i,j) = Hlos_fromEachLed(i,j,3);
        Hlos_led4(i,j) = Hlos_fromEachLed(i,j,4);
        %   method 1
%         channel(:,1) = conv(led_filter,Hlos_led1(i,j)*Ts*G_Con);
%         channel(:,2) = conv(led_filter,Hlos_led2(i,j)*Ts*G_Con);
%         channel(:,3) = conv(led_filter,Hlos_led3(i,j)*Ts*G_Con);
%         channel(:,4) = conv(led_filter,Hlos_led4(i,j)*Ts*G_Con);
%         
%         rxSignal(1,:) = filter(1,channel(:,1),txSignal(1,:));
%         rxSignal(2,:) = filter(1,channel(:,2),txSignal(2,:));
%         rxSignal(3,:) = filter(1,channel(:,3),txSignal(3,:));
%         rxSignal(4,:) = filter(1,channel(:,4),txSignal(4,:));
        % =================================================================
        %   method 2
        rxSignal(1,:) = conv(led_filter,txSignal(1,:));
        rxSignal(2,:) = conv(led_filter,txSignal(2,:));
        rxSignal(3,:) = conv(led_filter,txSignal(3,:));
        rxSignal(4,:) = conv(led_filter,txSignal(4,:));
        % =================================================================
        
        rx = sum(rxSignal);
        %rx = rx.*Ts.*G_Con;
        rx = awgn(rx, SNR, 'measured');
        %figure(2);
        [pxxRx, freq] = pwelch(rx,hanning(nfft),nfft/2,nfft,Fs);
        pxxRx_temps(i,j,:) = pxxRx;
        %plot(freq,smooth(10*log10(pxx2)),'*-','linewidth',1.2);
        %Pr(i,j,:) = [pxxRx(7),pxxRx(20),pxxRx(33),pxxRx(46)];
        Pr(i,j,:) = [pxxRx(14),pxxRx(40),pxxRx(66),pxxRx(96)];
        %   method 1
%         Pr_led1(i,j) = Pr(i,j,1);
%         Pr_led2(i,j) = Pr(i,j,2);
%         Pr_led3(i,j) = Pr(i,j,3);
%         Pr_led4(i,j) = Pr(i,j,4);
        
        % =================================================================
        %   method 2
        Pr_led1(i,j) = Pr(i,j,1)*Hlos_led1(i,j)*Ts*G_Con;
        Pr_led2(i,j) = Pr(i,j,2)*Hlos_led2(i,j)*Ts*G_Con;
        Pr_led3(i,j) = Pr(i,j,3)*Hlos_led3(i,j)*Ts*G_Con;
        Pr_led4(i,j) = Pr(i,j,4)*Hlos_led4(i,j)*Ts*G_Con;
        
               
        [de_a, de_b, de_c, de_d] = D_estimated(Adet,m,Pr_led1(i,j), Pr_led2(i,j), Pr_led3(i,j),Pr_led4(i,j), [Pt(1),Pt(2),Pt(3),Pt(4)],(room.height-receiver_height));
        re_a = sqrt(de_a^2-(room.height-receiver_height)^2);
        re_b = sqrt(de_b^2-(room.height-receiver_height)^2);
        re_c = sqrt(de_c^2-(room.height-receiver_height)^2);
        re_d = sqrt(de_d^2-(room.height-receiver_height)^2);
        
        A = [led.position_x(2)-led.position_x(1),led.position_y(2)-led.position_y(1);...
            led.position_x(3)-led.position_x(1),led.position_y(3)-led.position_y(1);...
            led.position_x(4)-led.position_x(1),led.position_y(4)-led.position_y(1)];
        B = [1/2*((re_a)^2-(re_b)^2) + 1/2*(led.position_x(2)^2+led.position_y(2)^2) - 1/2*(led.position_x(1)^2+led.position_y(1)^2);...
            1/2*((re_a)^2-(re_c)^2) + 1/2*(led.position_x(3)^2+led.position_y(3)^2) - 1/2*(led.position_x(1)^2+led.position_y(1)^2);...
            1/2*((re_a)^2-(re_d)^2) + 1/2*(led.position_x(4)^2+led.position_y(4)^2) - 1/2*(led.position_x(1)^2+led.position_y(1)^2)];
        xe = pinv(A' * A) * A' * B;
        
        locx(i,j) = xe(1);
        locy(i,j) = xe(2);
        locz(i,j) = receiver_height;
        
        error(i,j) = sqrt((X(i,j) - locx(i,j))^2 + (Y(i,j) - locy(i,j))^2 + (receiver_height - locz(i,j))^2);
    end
end



array = Pr_led1;
x = 0;
y = 1;
% Normalize to [0, 1]:
m = min(min(array));
range = max(max(array)) - m;
array = (array - m) ./ range;

% Then scale to [x,y]:
range2 = y - x;
normalized = (array.*range2) + x;


figure('Color','w');
s = surf(X,Y,normalized);%mechc xy坐标平面上有条纹
s.FaceColor = 'flat';
title('The normalized RSS values for all points on the receiver plan')
colorbar;
%caxis([0 0.4]);
xlabel('X(m)');
ylabel('Y(m)');
%view(2);



figure('Color','w');
%set(gca,'Color','w');
s = mesh(X,Y,error);%mechc xy坐标平面上有条纹
s.FaceColor = 'flat';
title('error')
colorbar;
%caxis([0 0.4]);
xlabel('X(m)');
ylabel('Y(m)');
%view(2);

figure('Color','w')
plot(X,Y,'bo',locx,locy,'*r','LineWidth',1.5);

figure('color','w');
stem3(X, Y, error,'filled','Color','k','LineWidth',1.5)
title('PE values at different points for the SNR of 20 dB')
xlabel('X(m)');
ylabel('Y(m)');
zlabel('PE(m)');

figure('color','w');
h0 = cdfplot(error(:).*100);%hold on
%h1 = cdfplot(error(2,:)*100);hold off
set(h0, 'color', 'g', 'LineStyle', '-.', 'LineWidth', 2);
%set(h1, 'color', 'b', 'LineStyle', '-.', 'LineWidth', 2);
%legend('SNR = 50dB','SNR = 55dB')
xlabel('Positioning Error (cm)');
ylabel('The CDF of the Positioning Error');

fprintf('SNR = %5.2f dB\n', SNR);
fprintf('Max error = %5.2f m\n', max(error(:)));
fprintf('Min error = %5.2f m\n', min(error(:)));
fprintf('Mean error = %5.2f m\n', mean(error(:)));
        

X_ind = 5;
Y_ind = 13;
pxxRx_led1 = pxxRx_temps(X_ind,Y_ind,:).*Hlos_led1(X_ind,Y_ind);
pxxRx_led2 = pxxRx_temps(X_ind,Y_ind,:).*Hlos_led2(X_ind,Y_ind);
pxxRx_led3 = pxxRx_temps(X_ind,Y_ind,:).*Hlos_led3(X_ind,Y_ind);
pxxRx_led4 = pxxRx_temps(X_ind,Y_ind,:).*Hlos_led4(X_ind,Y_ind);

pxxRx_leds = [pxxRx_led1(1,1:26),pxxRx_led2(1,27:52),pxxRx_led3(1,53:77),pxxRx_led4(1,78:end)];

figure
set(gca,'Color','w');
plot(freq,smooth(10*log10(pxxRx_leds)),'linewidth',1.2);
grid on
title('PSD of the received signal with E_b/N_0 = 20dB, receiver vertically below LED4')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
