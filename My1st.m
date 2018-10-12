%--------------------------------------------------------------------------
%This Program is the simulation of a QPSK communication system if the 
%channel is AWGN. We draw charachteristic curve of bit error rate(BER) of 
%the system versus Eb/No and compare with the BER of the theoritical data.

clc
clear all
close all
NumBits=1e6;
%Generating random data 
DataRn=round(rand(1,NumBits));
DataIn=sign(DataRn);
EbN0=(0:10);
BER=zeros(1,length(EbN0));  

%---Transmitter---------
%Gray mapping of bits into symbols
Column=length(DataIn)/2;
SigI=zeros(1,Column);
SigQ=SigI;

%using the formula for mapping to qpsk
QPSKMod=(-2)*sign(DataIn)+1;    
SigI=QPSKMod(1:2:length(DataIn)-1);
SigQ=QPSKMod(2:2:length(DataIn));
             
% complex mapped data
DataCx=SigI+1i*SigQ; 

%constellation diagram for QPSK
    %         ^
%         |
%  10 x   |   x 00   (odd bit, even bit)
%         |
%  -------+------->  real part (I channel)
%         |
%  11 x   |   x 01
%         |
% constellation diagram of the mapped data
%scatterplot(DataCx); 
%title('Signal diagram after qpsk mapped without noise');
%hold on ;
%grid on;

L=1;
while L <=length(EbN0)
    %--------CHANNEL-----------
    %Random noise generation and addition to the signal
    SNR=(10).^((EbN0(L))/10);
    Var=1/sqrt(2.*SNR);
    Noise=Var*(randn(1,length(DataCx))+ 1i*randn(1,length(DataCx))); 
    DataNo=DataCx+Noise;
    
    %--------RECEIVER----------
    %Demapped data and error calculations
    RX_Bits=zeros(1,NumBits);
   
    %for x=-2*I+1; where the inverse formula is I=-(x-1)/2
    D=1;
    for n=1:NumBits/2
       RX_Bits(D)= -(sign(real( DataNo(n)))-1)/2;
       RX_Bits(D+1)= -(sign(imag( DataNo(n)))-1)/2;
       D=D+2;
    end
    %Error Calculation
     NumbErrors=length(find(DataIn-RX_Bits));
     BER(L)=NumbErrors/NumBits;
     L=L+1;
end 

% Theoritical Bit error rate
TheoBER = 0.5*erfc(sqrt(10.^(EbN0/10)));

figure
semilogy(EbN0,TheoBER,'rs--','LineWidth',1.5,...
                'MarkerEdgeColor','k',...
                'MarkerSize',5);
hold on
semilogy(EbN0,BER,'gx-','LineWidth',1.5,...
                'MarkerEdgeColor','b',...
                'MarkerSize',5);
grid on;
axis([0 10 10^-5 1]);
xlabel('Eb/No,(dB)')
ylabel('Bit Error Rate (BER)')
title('Uncoded BER performance for QPSK in AWGN')
legend('Theoretical','Simulation')