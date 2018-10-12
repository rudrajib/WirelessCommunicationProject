%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Simulation of Communication System usng QPSK modulation under Rayleigh Fading
%                                                                                            
%                               Punit Godhani-  9709894   
%                               Balal Asar-     9755535
%                               Bilal Akhtar-   5208386
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
tic
N=1034;               % No of bits
EbN0_dB=1:25;         % Eb/N0 values in dB
E=1;                  % Normalization of Symbol Energy

% Loop for calculating Probability of error for Different values of Eb/No
for z=1:length(EbN0_dB)
    snr=10.^(EbN0_dB(z)/10); % Conversion of Eb/N0(dB) to linear value
    sgma=sqrt(E/(2*snr));    % Standard Deviation of Noise
    bit_error=0;             % Initializing Error counter
    bits=0;                  % Initializing Bit counter
    
    while bit_error<1e4
        data=round(rand(1,N)); % Generation of data
        
        data_I=data(1:2:end); %Inphase components of Data
        data_Q=data(2:2:end); %Quadrature components of Data
        
        % Modulation
        data_I_mod=data_I*(-2)+1;
        data_Q_mod=data_Q*(-2)+1;
        
        % Transmitted Signal
        Tx_sig=data_I_mod+sqrt(-1)*data_Q_mod;
        
        % Generation of AWGN Gaussian noise
        noise=sgma*(randn(1,517)+sqrt(-1)*randn(1,517));
        
        % Generation of Fading Coefficient
        h=sqrt(0.5)*randn+ sqrt(0.5)*sqrt(-1)*randn;
        
        %Received Signal under fading and AWGN noise
        Rx_sig= Tx_sig.*h+noise;
        
        % Dividing Received Signal by Fading Coefficient
        divide_h=Rx_sig./h;
        
        % Detection
        received_sig_I=real(divide_h)<=0;
        received_sig_Q=imag(divide_h)<=0;
        
        % Estimation of Data
        rx_sig= vertcat(received_sig_I,received_sig_Q);
        rx_msg=reshape(rx_sig,1,1034);
        
        bit_error=sum(xor(data,rx_msg))+bit_error; % Error Detection by comparing Estimated Data with the Transmitted Data
        bits=bits+N; %Incrementing the bit counter by number of bits
    end
    
    theoretical_ber(z)=1/2*(1-sqrt(snr/(1+snr))); % Calculation of Theoretical Bit Error Rate
    estimated_ber(z)=bit_error/bits; % Calculation of Estimated Bit Error Rate
end

% Plotting Bit Error Probability Vs Eb/N0 curve
close all
semilogy(EbN0_dB,theoretical_ber,'bs-','LineWidth',2);
hold on
semilogy(EbN0_dB,estimated_ber,'rx--','LineWidth',2);
grid on;
axis([1 max(EbN0_dB) 10^-5 1]); % Defining the Axis Range
legend('Theoretical', 'Simulation');
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Probability (Pb)');
title('Bit Error Probability (Pb) Vs Eb/N0 curve for QPSK with Rayleigh Fading');
toc
