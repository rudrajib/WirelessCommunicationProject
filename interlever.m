%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation of BCH (15,11,1) coded QPSK Communication System under Rayleigh Fading using Interleaver
%                                                                                            
%                      Punit Godhani-  9709894   
%                      Balal Asar-     9755535
%                      Bilal Akhtar-   5208386
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
tic
N=2326500;      % No of bits
EbN0_dB=0:25;   % Eb/N0 values in dB
E=1;            % Normalization of Symbol Energy
parmat=[1 1 1 1;0 1 1 1;1 0 1 1;1 1 0 1;1 1 1 0;0 0 1 1; 0 1 0 1;0 1 1 0;1 0 1 0; 1 0 0 1;1 1 0 0]; % Parity Matrix
Genmat=[parmat,eye(11)]; % Generator Matrix
H=[eye(4),parmat'];      % Parity Check Matrix
err=vertcat(zeros(1,15),eye(15));   % Error Pattern Table
S=err*H';     % Syndrome Table

% Loop for calculating Probability of error for Different values of Eb/No
for z=1:length(EbN0_dB)
    snr=10^(EbN0_dB(z)/10); % Conversion of Eb/N0(dB) to linear value
    sgma=sqrt((E/(2*snr))*(15/11));  % Standard Deviation of Noise
    bit_error=0;             % Initializing Error counter
    bits=0;                  % Initializing Bit counter
    x=1;
    y=1410*11;
    bits_gen=round(rand(1,N));   % Data Generation
    
    % Loop for calculating BER for each block in the Message
    % Depth of Interleaver = One Block of slow fading in BCH Coded system
    for i=1:N/(1410*11)
        data=bits_gen(x:y);
        word=reshape(data,11,1410).';   % Conversion of Data vector to 11 bit messages
        code_word=mod(word*Genmat,2);   % Generation of Code Word
        
        code_wordI=code_word(1:2:end,:); % Inphase Components of Code Word
        code_wordQ=code_word(2:2:end,:); % Quadrature Components of Code Word
        
        % Modulation of Code Word
        code_wordI_map=code_wordI*(-2)+1;
        code_wordQ_map=code_wordQ*(-2)+1;
        
        % Transmitted Signal
        Tx_code_word=code_wordI_map+sqrt(-1)*code_wordQ_map;
        
        % Generation of Fading Coefficients
        h=sqrt(0.5)*(randn(1,15)+sqrt(-1)*randn(1,15));
        
        % Generation of AWGN noise
        noise=sgma*(randn(705,15)+sqrt(-1)*randn(705,15));
        
        % Loop for Multiplying Fading Coefficients to Each bits in column
        for a=1:15
            rayleigh_sig(:,a)=h(a).*Tx_code_word(:,a);
        end
        
        % Received Signal Through AWGN channel
        Rx_code_word=rayleigh_sig+noise;
        
        % loop for Dividing each column by respective fading coefficient
        for a=1:15
            divide_h(:,a)=Rx_code_word(:,a)./h(a);
        end
        
        %Detection
        Rx_data_I=real(divide_h)<=0;
        Rx_data_Q=imag(divide_h)<=0;
        Rx_data_I_msg=reshape(Rx_data_I,1,705*15);
        Rx_data_Q_msg=reshape(Rx_data_Q,1,705*15);
        
        %Estimated Code Word
        Estimated_data=vertcat(Rx_data_I_msg,Rx_data_Q_msg);
        Estimate_data=reshape(Estimated_data,1410,15);
        
        % Syndrome Decoding
        s=mod(Estimate_data*H',2);
        
        % Error Detection and Correction
        for row=1:1410
            for counter=1:16
                if s(row,:)==S(counter,:)
                    corrected_word(row,:)=xor(Estimate_data(row,:),err(counter,:));
                end
            end
        end
        
        % Conversion of Code Word to Message Vector
        corrected_msg=corrected_word(:,5:15);
        corrected_msg_vector=reshape(corrected_msg',1,1410*11);
        
        bit_error=sum(xor(corrected_msg_vector,data))+bit_error; % Calculation of No. of Errors
        x=x+1410;
        y=y+1410;
    end
    % Calculation of Estimated Bit Error Rate using an interleaver
    estimated_ber(z)=bit_error/N;
    
    % Calculation of Theoretical Bit Error Rate for BCH coded under Rayleigh Fading
    theoretical_ber(z)=1/(4*snr*11/15);
    
end
% Plotting Bit Error Probability Vs Eb/N0 curve
close all
semilogy(EbN0_dB,theoretical_ber,'bs-','LineWidth',2);
hold on
semilogy(EbN0_dB,estimated_ber,'rx--','LineWidth',2);
grid on;
axis([0 max(EbN0_dB) 10^-5 1]); % Defining the Axis Range
legend('BER of BCH coded QPSK system', 'BER Using an Interleaver');
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Probability (Pb)');
title('Bit Error Probability (Pb) Vs Eb/N0 curve for BCH coded QPSK system with Interleaver under fading');
toc



