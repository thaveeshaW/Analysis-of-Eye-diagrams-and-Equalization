% Simulation Assignment – Eye diagrams and Equalization
% Wathudura T.R. - 210682D
% Dodangoda D.K.S.J. - 210150V

% Task 3
% Zero-forcing (ZF) equalizer for a 3-tap Multipath Channel

clear all; 
close all; 
clc;

% System Parameters
BitLength = 10^3;

% Generate a random binary sequence
BinarySequence = randi([0,1],1,BitLength);

% 2-PAM Signal
PAMsignal = 2 * BinarySequence - 1;
figure
stem(PAMsignal,'Marker','none');
title('2-PAM Signal');
grid


% Generate the received signal samples by convolving with channel response
h = [0.3 0.7 0.4];
ReceivedSignal = conv(PAMsignal,h);

BitErrors = zeros(10,5);
value = zeros(1,BitLength);

for n = 1:4
    % The matrix to obtain the ZF equalizer
    SignalWithNoise = awgn(ReceivedSignal,0);
    
    ShiftedMatrix = toeplitz([h([2:end]) zeros(1, 2 * n -1)],[h([2:-1:1]) zeros(1,2 * n - 1)]); 
    RequiredMatrix = zeros(1,2 * n + 1);
    RequiredMatrix(n + 1) = 1;
    ZF_equalizerCoefficients = [inv(ShiftedMatrix)*RequiredMatrix.'].';
    
    for k = 1:10
        % Add White Gaussian Noise
        SignalWithNoise = awgn(ReceivedSignal,k);

        % Convolve with filter to get response
        y_filtered = conv(SignalWithNoise,ZF_equalizerCoefficients);
        y_filtered = y_filtered(n+2:n+1+BitLength);% Consider the appropriate samples
            for i = 1:length(PAMsignal)
                if (y_filtered(i) > 0)
                    value(i) = 1;
                else
                    value(i) = 0;
                end
            end
        NoOfBitErrors = sum(value ~= BinarySequence);
        BitErrors(k,n) = NoOfBitErrors / BitLength;
    end
end

% AWGN Channel (Only noise)
for k = 1:10
    SignalWithNoise = awgn(PAMsignal,k);
     for i = 1:length(PAMsignal)
        if (SignalWithNoise(i) > 0) 
            value(i) = 1;
        else
            value(i) = 0;
        end
     end
     NoOfBitErrors = sum(value ~= BinarySequence);
     BitErrors(k,5) = NoOfBitErrors / BitLength;
end

% BER vs SNR graphs
figure
semilogy([1:10],BitErrors(:,1),'b*-','Linewidth',1.8);
hold on
semilogy([1:10],BitErrors(:,2),'g*-','Linewidth',1.8);
semilogy([1:10],BitErrors(:,3),'k*-','Linewidth',1.8);
semilogy([1:10],BitErrors(:,4),'m*-','Linewidth',1.8);
semilogy([1:10],BitErrors(:,5),'r*-','Linewidth',1.8);
axis([0 10 10^-3 0.5])
grid on
legend('3-Tap', '5-Tap','7-Tap','9-Tap','Noise Only');
xlabel('Eb/No in dB');
ylabel('BER');
title('BER vs Eb/No with ZF equalizer');