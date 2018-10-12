clc;
clear;
NumBits=1000000;              %  Number of bits 
DataIn=round(rand(1,NumBits));   %  To generate input bits of length NumBits, later must divide it to imaginary I and Real Q
signal=-1+2.*sign(DataIn);       %  Received =0.5*(1+sign(signal));
Min=0;
Max=10;
EbN0=(Min:Max);  %  Range of signal to Noise ratio, Eb/N0, horizontally
   
BER=zeros(1,length(EbN0));   %  Bit Error Rate
                             %  Transmitter side  
                             %  Gray mapping of bits into symbols
Column=length(DataIn)/2;
SigI=zeros(1,Column);
SigQ=SigI;
SigI=DataIn(1:2:NumBits-1);
SigQ=DataIn(2:2:NumBits);
SigI= -2.*SigI+1;
SigQ= -2.*SigQ+1;
symbole=SigI+1i.*SigQ;
L=1;
while L <=length(EbN0)                     
              % AWGN Addetive White Guassian Noise,to generate Random Noise and add to complex signal
    NoisS=10.^(EbN0(L)/10);
    NoisVariable=1/sqrt(2.*NoisS);
    RxSymbole=symbole+(NoisVariable*randn(1,length(symbole))+ 1i*NoisVariable*randn(1,length(symbole)));
                               %  At the Reciver side:
    RX_Bits=zeros(1,NumBits);
                              % To de-Map:
    D=1;
    for n=1:NumBits/2
       RX_Bits(D)= -(sign(real( RxSymbole(n)))-1)/2;
       RX_Bits(D+1)= -(sign(imag( RxSymbole(n)))-1)/2;
       D=D+2;
     end
                            % To Simulate Bit Error Rate:
     NumbErrors=find(xor(RX_Bits,DataIn));
     NumbErrors=size(NumbErrors,2);
     BER(L)=NumbErrors/NumBits;
     L=L+1;
end     %  end of while loop
                          % To plot and lable on the semilog graphs
Theory_BER=0.5*erfc(sqrt(10.^(EbN0/10)));  % Theroretical Values
semilogy(EbN0,Theory_BER,'ro-');  % Theoretical graph, Uncoded
hold on
semilogy(EbN0,BER,'bx-');         % Simulation graph
xlabel('Eb/No,(dB)')
ylabel('Bit Error Rate')
title('QPSK Performance, Uncoded')
legend('Theoretical','Simulation')
grid on