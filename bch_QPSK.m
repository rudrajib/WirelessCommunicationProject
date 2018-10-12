%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Simulation of BCH(15,11,1)coded Communication System using QPSK modulation
%                                                                          
%                               Punit Godhani-  9709894   
%                               Balal Asar-     9755535
%                               Bilal Akhtar-   5208386
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
tic
N=99000;  % No of bits
EbN0_dB=0:10;   % Eb/N0 values in dB
E=1;      % Normalization of Symbol Energy
parmat=[1 1 1 1;0 1 1 1;1 0 1 1;1 1 0 1;1 1 1 0;0 0 1 1; 0 1 0 1;0 1 1 0;1 0 1 0; 1 0 0 1;1 1 0 0]; % Parity Matrix
Genmat=[parmat,eye(11)]; % Generator Matrix
H=[eye(4),parmat'];   % Parity Check Matrix
err=vertcat(zeros(1,15),eye(15));   % Error Pattern Table
S=err*H';      % Syndrome Table

% Loop for calculating Probability of error for Different values of Eb/No
for z=1:length(EbN0_dB)
    snr=10^(EbN0_dB(z)/10); % Conversion of Eb/N0(dB) to linear value
    sgma=sqrt((E/(2*snr))*(15/11));  % Standard Deviation of Noise
    bit_error=0;             % Initializing Error counter
    
    data=round(rand(1,N));   % Data Generation
    word=reshape(data,11,N/11).'; % Conversion of Data vector to 11 bit messages
    code_word=reshape(mod(word*Genmat,2).',1,(N*15/11));  % Generation of Code Word
    
    code_wordI=code_word(1:2:end); % Inphase Components of Code Word
    code_wordQ=code_word(2:2:end); % Quadrature Components of Code Word
    
    % Modulation of Code Word
    code_wordI_mod=code_wordI*(-2)+1;
    code_wordQ_mod=code_wordQ*(-2)+1;
    
    % Transmitted Signal
    Tx_code_word=code_wordI_mod+sqrt(-1)*code_wordQ_mod;
    
    % Generation of AWGN noise
    noise=sgma*(randn(1,length(code_word)/2)+sqrt(-1)*randn(1,length(code_word)/2));
    
    % Received Signal Through AWGN channel
    Rx_code_word=Tx_code_word+noise;
    
    %Detection
    Rx_data_I=real(Rx_code_word)<=0;
    Rx_data_Q=imag(Rx_code_word)<=0;
    
    %Estimated Code Word
    Estimated_data=vertcat(Rx_data_I,Rx_data_Q);
    estimate_code_word=reshape(Estimated_data,15,length(code_word)/15).';
    
    % Syndrome Decoding
    s=mod(estimate_code_word*H',2);
    
    % Error Detection and Correction
    for row=1:length(code_word)/15
        for counter=1:16
            if s(row,:)==S(counter,:)
                corrected_word(row,:)=xor(estimate_code_word(row,:),err(counter,:));
            end
        end
    end
    
    % Conversion of Code Word to Message Vector
    corrected_msg=corrected_word(:,5:15);
    corrected_vector=reshape(corrected_msg',1,N);
    
    % Calculation of Simulated Bit Error Rate
    bit_error=sum(xor(corrected_vector,data));
    estimated_ber(z)=bit_error/(N);
    
    % Calculation of Theoretical Bit Error Rate
    pr_code=(1/2)*erfc(sqrt(snr*(11/15)));
    funct=0;
    for i=2:15
        funct=i*nchoosek(15,i)*pr_code^i*(1-pr_code)^(15-i)+funct;
        theoretical_ber(z)=(1/15)*funct;
    end
    
end

% Plotting Bit Error Probability Vs Eb/N0 curve
close all
semilogy(EbN0_dB,theoretical_ber,'bs-','LineWidth',2);
hold on
semilogy(EbN0_dB,estimated_ber,'rx--','LineWidth',2);
grid on;
axis([0 max(EbN0_dB) 10^-5 1]); % Defining the Axis Range
legend('Theoretical', 'Simulation');
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Probability (Pb)');
title('Bit Error Probability (Pb) Vs Eb/N0 curve for BCH coded QPSK system');
toc



