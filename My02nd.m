%--------------------------------------------------------------------------
%This Program is the simulation of a QPSK communication system if the 
%channel is AWGN. We consider single error correcting (15,11) code and 
%draw charachteristic curves of bit error rate(BER) of the system versus 
%Eb/No and compare with the BER of the uncoded system

clc
clear all
close all
NumBits=11e5; %matching bit rate

%Generating random data 
DataRn=round(rand(1,NumBits));
DataIn=sign(DataRn);
EbN0=(0:10);  
BER=zeros(1,length(EbN0));  
BER2=zeros(1,length(EbN0)); 

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
S=mod(E*Ht,2);          % Syndrome table Se
CWord=zeros(1,length(DataIn)*n/k);

% Conversion 11 bit messages to code word
CWord=reshape(DataIn,11,NumBits/11).'; 
CodedData=reshape(mod(CWord*G,2).',1,(NumBits*15/11));                    

%------Transmitter side--------                              
%Gray mapping of bits into symbols
Column=length(CodedData)/2;
SigI=zeros(1,Column);
SigQ=SigI;

%Gray mapping of bits into symbols
Col=length(DataIn)/2;
I=zeros(1,Col);
Q=I;

%using the formula for mapping to qpsk
CodedData=(-2)*sign(CodedData)+1;
SigI=CodedData(1:2:length(CodedData)-1);
SigQ=CodedData(2:2:length(CodedData));

%using the formula for mapping to qpsk
QPSKDMod=(-2)*sign(DataIn)+1;    
I=QPSKDMod(1:2:length(DataIn)-1);
Q=QPSKDMod(2:2:length(DataIn));

% complex mapped data
DataCx=SigI+1i.*SigQ;
DataUnCx=I+1i*Q;

LL=1;
while LL <= length(EbN0) 

    SNR2=(10).^((EbN0(LL))/10);
    Var2=sqrt(1/(2*SNR2));
    Noise2=Var2*(randn(1,length(DataUnCx))+ 1i*randn(1,length(DataUnCx))); 
    DataUnNo=DataUnCx+Noise2;
    
     %--------RECEIVER----------
    %Demapped data and error calculations
    RX_Bits2=zeros(1,length(DataIn));
   
    %for x=-2*I+1; where the inverse formula is I=-(x-1)/2
    R=1;
    for J=1:NumBits/2
       RX_Bits2(R)= -(sign(real( DataUnNo(J)))-1)/2;
       RX_Bits2(R+1)= -(sign(imag( DataUnNo(J)))-1)/2;
       R=R+2;
    end
    %Error Calculation
     NumbErrors=length(find(DataIn-RX_Bits2));
     BER2(LL)=NumbErrors/NumBits;
     LL=LL+1;
end
               
semilogy(EbN0,BER2,'rx-','LineWidth',1.5,'MarkerEdgeColor','k',...
                  'MarkerSize',5);  
hold on  
L=1;
while L <= length(EbN0)                   
              
    %--------CHANNEL-----------
    %Random noise generation and addition to the signal
    SNR=10.^(EbN0(L)/10);
    Var=sqrt((1/(2*SNR))*(15/11));  % Standard Deviation of Noise
    Noise=Var*(randn(1,length(DataCx))+ 1i*randn(1,length(DataCx))); 
    DataNo=DataCx+Noise;
    
     
    %--------RECEIVER----------
    %Demapped data and error calculations
    RX_Bits=zeros(1,length(CodedData));
    
    x=1;
    for func=1:length(DataNo)
        RX_Bits(x)= -(sign(real(DataNo(func)))-1)/2;
        RX_Bits(x+1)= -(sign(imag(DataNo(func)))-1)/2;
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

semilogy(EbN0,BER,'gx-','LineWidth',1.5,'MarkerEdgeColor','b',...
                  'MarkerSize',5);    
     
grid on   
axis([0 10 10^-5 1]);
xlabel('Eb/No,(dB)')
ylabel('Bit Error Rate (BER)')
title('Uncoded and Coded BER performance')
legend('Simulation-Uncoded','Simulation-Coded')