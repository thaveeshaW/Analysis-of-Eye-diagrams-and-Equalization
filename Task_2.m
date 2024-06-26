% Simulation Assignment – Eye diagrams and Equalization
% Wathudura T.R. - 210682D
% Dodangoda D.K.S.J. - 210150V

% Task 2
% With Additive White Gaussian Noise (AWGN)

clear all; 
close all; 
clc;

% System Parameters
BitLength = 10^3; % No of Bits Transmitted
SampleFreq = 10; % Sampling frequency (Hz)
time = -SampleFreq:1/SampleFreq:SampleFreq; % Time Array
SNR_dB = 10;
NoisePower = 1./(10.^(0.1*SNR_dB)); % Noise Power (Eb = 1 in BPSK)

% Generate the BPSK Signal (With noise)
% Map 0 -> -1 and 1 -> 1 (0-phase and 180-phase)
BPSKSignal = 2*(rand(1,BitLength)>0.5)-1;
t = 0:1/SampleFreq:99/SampleFreq;

% Noise Array Generation using SNR = 10dB
Noise1D = normrnd (0 , sqrt(NoisePower/2), [1, BitLength]);
AWGN_TX = BPSKSignal + Noise1D;
figure;
stem(t, AWGN_TX(1:100)); xlabel('Time'); ylabel('Amplitude');
title('BPSK Impulse Train (With Noise)');
axis([0 10 -1.5 1.5]); grid on;

% Sinc Pulse
Sinc_Num = sin(pi*time); % numerator sinc
Sinc_Den = (pi*time); % denominator sinc
Sinc_DenZero = find(abs(Sinc_Den) < 10^-10); % Find t=0 position
Sinc_Filt = Sinc_Num./Sinc_Den;
Sinc_Filt(Sinc_DenZero) = 1; % Define t=0 value

% Raised Cosine Pulse for 0.5 roll-off
roll_off = 0.5;
cos_Num = cos(roll_off*pi*time);
cos_Den = (1 - (2 * roll_off * time).^2);
cos_DenZero = abs(cos_Den)<10^-10;
RaisedCosine = cos_Num./cos_Den;
RaisedCosine(cos_DenZero) = pi/4;
RC_gamma5 = Sinc_Filt.*RaisedCosine; % Getting the complete raised cosine pulse

% Raised Cosine Pulse for 1 roll-off
roll_off = 1;
cos_Num = cos(roll_off * pi * time);
cos_Den = (1-(2 * roll_off * time).^2);
cos_DenZero = find(abs(cos_Den)<10^-20);
RaisedCosine = cos_Num./cos_Den;
RaisedCosine(cos_DenZero) = pi/4;
RC_gamma1 = Sinc_Filt.*RaisedCosine; % Get the complete raised cosine pulse

% Upsample the transmit sequence (With Noise)
AWGNTx_Upsample = [AWGN_TX;zeros(SampleFreq-1,length(BPSKSignal))];
AWGNTx_U = AWGNTx_Upsample(:);
figure;
stem(t, AWGNTx_U(1:100)); xlabel('Time'); ylabel('Amplitude');
title('Upsampled BPSK Impulse Train (With Noise)');
axis([0 10 -1.5 1.5]); grid on;

% Pulse Shaped sequences (With Noise)
Conv_sincnoise = conv(AWGNTx_U,Sinc_Filt);
Conv_RC5noise = conv(AWGNTx_U,RC_gamma5);
Conv_R1noise = conv(AWGNTx_U,RC_gamma1);

% Take only the first 10000 samples (With noise)
Conv_sincnoise = Conv_sincnoise(1:10000);
Conv_RC5noise = Conv_RC5noise(1:10000);
Conv_R1noise = Conv_R1noise(1:10000);

% Reshape the sequences to build Eye Diagrams (With Noise)
Conv_sincnoise_reshape = reshape(Conv_sincnoise, SampleFreq*2, BitLength*SampleFreq/20).';
Conv_RC5noise_reshape = reshape(Conv_RC5noise,SampleFreq*2,BitLength*SampleFreq/20).';
Conv_R1noise_reshape = reshape(Conv_R1noise,SampleFreq*2,BitLength*SampleFreq/20).';

% Plot the Eye Diagrams (With Noise)
% Eye diagram for Sinc pulse
figure;
plot(0:1/SampleFreq:1.99, Conv_sincnoise_reshape.', 'b');
title('Eye diagram for Sinc pulse (With Noise)');
xlabel('Time'); ylabel('Amplitude');
axis([0 2 -2.4 2.2]);
grid on

% Eye diagram for RC pulse with Gamma = 0.5
figure;
plot(0:1/SampleFreq:1.99, Conv_RC5noise_reshape.', 'b');
title('Eye diagram for RC pulse with Gamma = 0.5 (With Noise)');
xlabel('Time'); ylabel('Amplitude');
axis([0 2 -2.4 2.2]);
grid on

% Eye diagram for RC pulse with Gamma = 1
figure;
plot(0:1/SampleFreq:1.99, Conv_R1noise_reshape.', 'b');
title('Eye diagram for RC pulse with Gamma = 1 (With Noise)');
xlabel('Time'); ylabel('Amplitude');
axis([0 2 -2.4 2.2]);
grid on