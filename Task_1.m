% Simulation Assignment – Eye diagrams and Equalization
% Wathudura T.R. - 210682D
% Dodangoda D.K.S.J. - 210150V

% Task 1

clear all; 
close all; 
clc;

% System Parameters
BitLength = 10^3; % No of Bits Transmitted
SampleFreq = 10; % Sampling frequency (Hz)
time = -SampleFreq:1/SampleFreq:SampleFreq; % Time Array

% Generate the BPSK Signal
% Map 0 -> -1 and 1 -> 1 (0-phase and 180-phase)
BPSKSignal = 2*(rand(1,BitLength)>0.5)-1;
t = 0:1/SampleFreq:99/SampleFreq;
stem(t, BPSKSignal(1:100)); xlabel('Time'); ylabel('Amplitude');
title('BPSK Impulse Train');
axis([0 10 -1.2 1.2]); grid on;

% Sinc Pulse
Sinc_Num = sin(pi*time); % numerator sinc
Sinc_Den = (pi*time); % denominator sinc
Sinc_DenZero = find(abs(Sinc_Den) < 10^-10); % Find t=0 position
Sinc_Filt = Sinc_Num./Sinc_Den;
Sinc_Filt(Sinc_DenZero) = 1; % Define t=0 value
figure;
plot(time, Sinc_Filt);
title('Sinc Pulse');
xlabel('Time'); ylabel('Amplitude');
axis([-SampleFreq SampleFreq -0.5 1.2]); grid on

% Raised Cosine Pulse for 0.5 roll-off
roll_off = 0.5;
cos_Num = cos(roll_off*pi*time);
cos_Den = (1 - (2 * roll_off * time).^2);
cos_DenZero = abs(cos_Den)<10^-10;
RaisedCosine = cos_Num./cos_Den;
RaisedCosine(cos_DenZero) = pi/4;
RC_gamma5 = Sinc_Filt.*RaisedCosine; % Getting the complete raised cosine pulse
figure;
plot(time, RC_gamma5);
title('Raised Cosine Pulse for Gamma = 0.5');
xlabel('Time'); ylabel('Amplitude');
axis([-SampleFreq SampleFreq -0.5 1.2]); grid on

% Raised Cosine Pulse for 1 roll-off
roll_off = 1;
cos_Num = cos(roll_off * pi * time);
cos_Den = (1-(2 * roll_off * time).^2);
cos_DenZero = find(abs(cos_Den)<10^-20);
RaisedCosine = cos_Num./cos_Den;
RaisedCosine(cos_DenZero) = pi/4;
RC_gamma1 = Sinc_Filt.*RaisedCosine; % Get the complete raised cosine pulse
figure;
plot(time, RC_gamma1);
title('Raised Cosine Pulse for Gamma = 1');
xlabel('Time'); ylabel('Amplitude');
axis([-SampleFreq SampleFreq -0.5 1.2]); grid on

% Upsample the transmit sequence
BPSK_Upsample = [BPSKSignal;zeros(SampleFreq-1,length(BPSKSignal))]; % Upsample the BPSK to match the sampling frequency
BPSK_U = BPSK_Upsample(:).';
figure;
stem(t, BPSK_U(1:100)); xlabel('Time'); ylabel('Amplitude');
title('BPSK Impulse Train after upsampling');
axis([0 10 -1.2 1.2]); grid on;

% Pulse Shaped sequences
Conv_sincpulse = conv(BPSK_U, Sinc_Filt);
Conv_RCgamma5 = conv(BPSK_U,RC_gamma5);
Conv_RCgamma1 = conv(BPSK_U,RC_gamma1);


% Take only the first 10000 samples
Conv_sincpulse = Conv_sincpulse(1:10000);
Conv_RCgamma5 = Conv_RCgamma5(1:10000);
Conv_RCgamma1 = Conv_RCgamma1(1:10000);

% Reshape the sequences to build Eye Diagrams
Conv_sincpulse_reshape = reshape(Conv_sincpulse, SampleFreq*2, BitLength*SampleFreq/20).';
Conv_RCgamma5_reshape = reshape(Conv_RCgamma5,SampleFreq*2,BitLength*SampleFreq/20).';
Conv_RCgamma1_reshape = reshape(Conv_RCgamma1,SampleFreq*2,BitLength*SampleFreq/20).';

% Plot the Eye Diagrams
% Eye diagram for Sinc pulse
figure;
plot(0:1/SampleFreq:1.99, real(Conv_sincpulse_reshape).', 'b');
title('Eye diagram for Sinc pulse');
xlabel('Time'); ylabel('Amplitude');
axis([0 2 -2.4 2.2]);
grid on

% Eye diagram for RC pulse with Gamma = 0.5
figure;
plot(0:1/SampleFreq:1.99, Conv_RCgamma5_reshape.','b');
title('Eye diagram for RC pulse with Gamma = 0.5');
xlabel('Time'); ylabel('Amplitude');
axis([0 2 -2.5 2.5]);
grid on

% Eye diagram for RC pulse with Gamma = 1
figure;
plot(0:1/SampleFreq:1.99, Conv_RCgamma1_reshape.','b');
title('Eye diagram for RC pulse with Gamma = 1');
xlabel('Time'); ylabel('Amplitude');
axis([0 2 -1.5 1.5 ]);
grid on