%--------------------------------------------------------------------------
%In this part of the project we have simulate the communication system 
%having Doppler Spread due to mobile movement in the channel.we Simulate
%and draw the bit error rate of the system for receiver speed 60 Km/Hour
%--------------------------------------------------------------------------

clc
clear all
close all
NumBits=11e5;

%Generating random data 
DataRn=round(rand(1,NumBits));
DataIn=sign(DataRn);

%---Transmitter---------
%Gray mapping of bits into symbols
Column=length(DataIn)/2;
SigI=zeros(1,Column);
SigQ=SigI;

%using the formula for mapping to qpsk
QPSKData=(-2)*sign(DataIn)+1;    
SigI=QPSKData(1:2:NumBits-1);
SigQ=QPSKData(2:2:NumBits);
             
% complex mapped data
DataCx=SigI+1i*SigQ; 

c=3e8; %velocity of light
fc=10e9; %carrier frequency
lembda=c/fc;
v=(60*1000)/3600; %receiver velocity
fm=v/lembda;
Tc=9/(16*pi*fm); %Coherence Time
Tb=2/NumBits;   %Symbol Time
B=Tc/Tb;    %block
B=round(B);

%noise of fading channel
Ch=zeros(1,length(DataCx)+length(B));

for n=1:length(B):length(DataCx)
        
        Ch(n:n+length(B)-1)=0.5*(randn*ones(1,length(B))+1i*randn*ones(1,length(B)));
    end
ChMax=Ch(1:length(DataCx));   

L=1;
EbNo=(0:25);
while L <=length(EbNo)
    %--------CHANNEL-----------
    %Random noise generation and addition to the signal
    SNR=(10)^((EbNo(L))/10);
    Var=1/sqrt(2*SNR);
    Noise=Var*(randn(1,length(DataCx))+ 1i*randn(1,length(DataCx))); 
    DataNo=DataCx.*ChMax+Noise;
    
    %--------RECEIVER----------
    %Demapped data and error calculations
    DataNo=DataNo./ChMax;
    RX_Bits=zeros(1,NumBits);
   
    %for x=-2*I+1; where the inverse formula is I=-(x-1)/2
    D=1;
    for n=1:NumBits/2
       RX_Bits(D)= -(sign(real( DataNo(n)))-1)/2;
       RX_Bits(D+1)= -(sign(imag( DataNo(n)))-1)/2;
       D=D+2;
    end
    
     NumbErrors=length(find(DataIn-RX_Bits));
     BER(L)=NumbErrors/NumBits;
     L=L+1;
end 

% Theoritical Bit error rate
ebn=10.^(EbNo/10);
TheoBER=.5.*(1-sqrt(ebn./(1+ebn)));


figure
semilogy(EbNo,TheoBER,'rs--','LineWidth',1.5,...
                'MarkerEdgeColor','k',...
                'MarkerSize',5);
hold on
semilogy(EbNo,BER,'gx-','LineWidth',1.5,...
                'MarkerEdgeColor','b',...
                'MarkerSize',5);

grid on;
axis([0 25 10^-5 1]);
xlabel('Eb/No,(dB)')
ylabel('Bit Error Rate')
title('Uncoded BER for QPSK in Rayheigh Fading')
legend('Theoretical','Simulation','No fading')