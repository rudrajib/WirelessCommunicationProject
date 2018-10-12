%--------------------------------------------------------------------------
%This Program is the simulation of a QPSK communication system if the 
%channel is AWGN. We consider single error correcting (15,11) code and 
%draw charachteristic curves of bit error rate(BER) of the system versus 
%Eb/No and compare with the BER of the theoritical data.

clc
clear all
close all
NumBits=11e5; %matching bit rate

%Generating random data 
DataRn=round(rand(1,NumBits));
DataIn=sign(DataRn);
EbN0=(0:10);  
BER=zeros(1,length(EbN0));   

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

%using the formula for mapping to qpsk
CodedData=(-2)*sign(CodedData)+1;
SigI=CodedData(1:2:length(CodedData)-1);
SigQ=CodedData(2:2:length(CodedData));

% complex mapped data
DataCx=SigI+1i.*SigQ;

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

% Calculation of Theoretical Bit Error Rate
TheoBER =0.5*erfc(sqrt(10.^(EbN0/10)*(11/15)));
i=1;
t=1; %Error Correction Capablity
TheoPB=zeros(1,length(EbN0));
while i<=length(EbN0)
    for j=t+1:n
        TheoPB(i)=TheoPB(i)+j*nchoosek(n,j)*TheoBER(i).^j*(1-TheoBER(i)).^(n-j);
    end
    TheoPB(i)=TheoPB(i)/n;
    i=i+1;
end        

semilogy(EbN0,TheoPB,'rs--','LineWidth',1.5,'MarkerEdgeColor','k',...
                  'MarkerSize',5);   
hold on
semilogy(EbN0,BER,'gx-','LineWidth',1.5,'MarkerEdgeColor','b',...
                  'MarkerSize',5);                     
grid on   
axis([0 10 10^-5 1]);
xlabel('Eb/No,(dB)')
ylabel('Bit Error Rate (BER)')
title('Coded BER performance for QPSK in AWGN')
legend('Theoretical','Simulation')