%--------------------------------------------------------------------------
%This Program is the simulation of a QPSK communication system if the 
%channel is AWGN under rayheigh fading. We consider single error correcting 
%(15,11) code and draw charachteristic curves of bit error rate(BER) of the 
%system versus Eb/No and compare with the BER of the uncoded system.

clc;
clear all;
close all;
NumBits=11e5;              %matching bit rate

%Generating random data 
DataRn=round(rand(1,NumBits));
DataIn=sign(DataRn);
EbNo=0:25;  
BER=zeros(1,length(EbNo));
BER2=zeros(1,length(EbNo));

%Gray mapping of bits into symbols
Column2=length(DataIn)/2;
SigI2=zeros(1,Column2);
SigQ2=SigI2;

%using the formula for mapping to qpsk
QPSKData=(-2)*sign(DataIn)+1;    
SigI2=QPSKData(1:2:NumBits-1);
SigQ2=QPSKData(2:2:NumBits);
             
% complex mapped data
DataCx2=SigI2+1i*SigQ2; 

%---BCH Coder------------
% using single error correcting (15,11) code  
n=15;
k=11;
In=eye(k);     % Idendity Matrix In
               % Parity Matrix P 
P=[1 1 1 1;0 1 1 1;1 0 1 1;1 1 0 1;1 1 1 0 ;0 0 1 1;0 1 0 1;...
   0 1 1 0;1 0 1 0;1 0 0 1;1 1 0 0];
Pt=P.';         % Pt is transpose matrix of matrix P
G=[P In];
H=[eye(n-k) Pt];
Ht=H.';         %  Ht is transpose matrix of Matrix H
E=[zeros(1,n);diag(ones(1,n))];  % Error Matrix E
S=mod(E*Ht,2);     % Syndrome table Se
Se=[S E];          
CWord=zeros(1,length(DataIn)*n/k);

% Conversion 11 bit messages and code word
CWord=reshape(DataIn,11,NumBits/11).'; 
CodedData=reshape(mod(CWord*G,2).',1,(NumBits*15/11)); 

%------Transmitter side--------      
%Gray mapping of bits into symbols
Column=length(CodedData)/2;
SigI=zeros(1,Column);
SigQ=SigI;

%using the formula for mapping to qpsk
CodedData=(-2)*sign(CodedData)+1;
SigI=CodedData(1:2:length(CodedData)-1);
SigQ=CodedData(2:2:length(CodedData));

% complex mapped data
DataCx=SigI+1i.*SigQ;

c=3e8; %velocity of light
fc=10e9; %carrier frequency
lembda=c/fc;
v=(60*1000)/3600; %receiver velocity
fm=v/lembda;
Tc=9/(16*pi*fm); %Coherence Time
Ts=1/NumBits;    %Symbol Time
Tb=Ts.*15/11;
B2=Tc/Ts;
B=Tc/Tb; %block
B2=round(B2);
B=round(B);

%noise of fading channel
Chh=zeros(1,length(DataCx2)+B2);

for nn=1:B2:length(DataCx2)
        
        Chh(nn:nn+B2-1)=0.5*(randn*ones(1,B2)+1i*randn*ones(1,B2));
    end
ChMax2=Chh(1:length(DataCx2));   

L2=1;
EbNo=(0:25);
while L2 <=length(EbNo)
    %--------CHANNEL-----------
    %Random noise generation and addition to the signal
    SNRR=(10)^((EbNo(L2))/10);
    Varr=1/sqrt(2*SNRR);
    Noisee=Varr*(randn(1,length(DataCx2))+ 1i*randn(1,length(DataCx2))); 
    DataNoo=DataCx2.*ChMax2+Noisee;
    
    %--------RECEIVER----------
    %Demapped data and error calculations
    DataNoo=DataNoo./ChMax2;
    RX_Bitss=zeros(1,NumBits);
   
    %for x=-2*I+1; where the inverse formula is I=-(x-1)/2
    DD=1;
    for nnn=1:NumBits/2
       RX_Bitss(DD)= -(sign(real( DataNoo(nnn)))-1)/2;
       RX_Bitss(DD+1)= -(sign(imag( DataNoo(nnn)))-1)/2;
       DD=DD+2;
    end
    
     NumbErrorss=length(find(DataIn-RX_Bitss));
     BER2(L2)=NumbErrorss/NumBits;
     L2=L2+1;
end 

semilogy(EbNo,BER2,'bx-','LineWidth',1.5,'MarkerEdgeColor','r',...
                 'MarkerSize',5);   
hold on

Ch=zeros(1,length(DataCx)+B);   
    for n1=1:B:length(DataCx)      
        Ch(n1:n1+B-1)=0.5*randn*ones(1,B)+j*0.5*randn*ones(1,B);
    end
    ChMax=Ch(1:length(DataCx)); 
    
L=1;    
while L <= length(EbNo)                   
    %--------CHANNEL-----------
    %Random noise generation and addition to the signal         
    SNR=10.^(EbNo(L)/10);
    Var=sqrt((11/15)/(2.*SNR));
    Noise=(Var*randn(1,length(DataCx))+ 1i*Var*randn(1,length(DataCx))); 
    DataNo=DataCx.*ChMax+Noise ;
    
    %--------RECEIVER----------
    %Demapped data and error calculations
    DataNo=DataNo./ChMax;
    RX_Bits=zeros(1,length(CodedData));
                                 % To Decode
    x=1;
    for j=1:length(DataNo)
        RX_Bits(x)= -(sign(real(DataNo(j)))-1)/2;
        RX_Bits(x+1)= -(sign(imag(DataNo(j)))-1)/2;
        x=x+2;
    end
    %Estimated Code Word
    RX_CodedBits=reshape(RX_Bits,15,length(CodedData)/15).';
    
    % Syndrome Decoding
    s=mod(RX_CodedBits*Ht,2);
    
    % Error Detection and Correction
    for Row=1:length(CodedData)/15
        for ConT=1:16
            if s(Row,:)==S(ConT,:)
                CorrBits(Row,:)=xor(RX_CodedBits(Row,:),E(ConT,:));
            end
        end
    end
    
    % Conversion of Message Vector from Code Word
    CorrMsgs=CorrBits(:,5:15);
    CorrVect=reshape(CorrMsgs',1,NumBits);
    
    % Simulated Bit Error Rate Calculation
    Errors=find(xor(CorrVect,DataIn));    
    Errors=size(Errors,2);
    BER(L)=Errors/NumBits;
    L=L+1; 
end   
% Calculation of Theoretical Bit Error Rate                             
Nois=10.^(EbNo/10);
TheoBER=0.5.*(1-sqrt(Nois./(1+Nois)));

semilogy(EbNo,BER,'gx-','LineWidth',1.5,'MarkerEdgeColor','r',...
                  'MarkerSize',5); 
hold on

semilogy(EbNo,TheoBER,'rs--','LineWidth',1.5,'MarkerEdgeColor','k',...
                  'MarkerSize',5);
grid on
axis([1 25 10^-5 1]);
xlabel('Eb/No,(dB)')
ylabel('Bit Error Rate')
title('Uncoded and Coded BER Under Rayleigh Fading')
legend('UnCoded','Coded','Theory')