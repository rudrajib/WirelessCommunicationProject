%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Simulation of Communication System using QPSK modulation
%                                                                          
%                               Punit Godhani-  9709894   
%                               Balal Asar-     9755535
%                               Bilal Akhtar-   5208386
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
tic
N=10^6;                 % No of bits
EbN0_dB=0:15;           % Eb/N0 values in dB
E=1;                    % Normalization of bit Energy

% Loop for calculating Probability of error for Different values of Eb/No
for z=1:length(EbN0_dB)
    snr=10^(EbN0_dB(z)/10); % Conversion of Eb/N0(dB) to linear value
    sgma=sqrt(E/(2*snr));   % Standard Deviation of Noise
    errors=0;  % Initializing Error counter
    
    % Bit Generation with equal probability
    Tx_bits=round(rand(1,N));
    Data_I=Tx_bits(1:2:end); % Inphase components of Data
    Data_Q=Tx_bits(2:2:end); %Quadrature Components of Data
    
    % Modulation of Data
    Data_I_map=Data_I*(-2)+1;
    Data_Q_map=Data_Q*(-2)+1;
    
    % Transmitted Signal
    Tx_sig=Data_I_map+sqrt(-1)*Data_Q_map;
    
    % Generation Of AWGN Noise
    noise=sgma*(randn(1,N/2)+sqrt(-1)*randn(1,N/2)); % White Gaussian noise with variance N0/2 per dimension
    
    % Received signal Through AWGN channel
    Rx_sig= Tx_sig+noise;
    
    %Detection
    Rx_data_I=real(Rx_sig)<=0;
    Rx_data_Q=imag(Rx_sig)<=0;
    
    % Estimated Data
    Estimated_data=vertcat(Rx_data_I, Rx_data_Q);
    Estimated_bits=reshape(Estimated_data,1,N);
    
    % Error Detection by comparing estimated data with the transmitted data
    errors=sum(xor(Estimated_bits,Tx_bits));
    
    % Counting the total number of bit errors by adding bit difference in each estimated symbol
    theoretical_ber(z)=(1/2)*erfc(sqrt(snr)); % Calculation of Theoretical Bit Error Rate
    estimated_ber(z)=errors/N; % Calculation of Estimated Bit Error Rate
end

% Plotting Bit Error Probability Vs Eb/N0 curve
semilogy(EbN0_dB,theoretical_ber,'bs-','LineWidth',2);
hold on
semilogy(EbN0_dB,estimated_ber,'rx--','LineWidth',2);
grid on;
axis([0 max(EbN0_dB) 10^-5 1]); % Defining the Axis Range
legend('Theoretical', 'Simulation');
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Probability (Pb)');
title('Bit Error Probability (Pb) Vs Eb/N0 curve for QPSK');
toc
