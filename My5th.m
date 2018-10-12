%--------------------------------------------------------------------------
%This Program is the simulation of a QPSK communication system if the 
%channel is AWGN under rayheigh fading. We consider single error correcting 
%(15,11) code and interleaver to improve the performance. We draw 
%charachteristic curves of BER of the system versus Eb/No and compare with 
%the BER of without interleaver and uncoded data.

clc;
clear all;
close all;
NumBits=242e3;   %matching bit rate

%Generating random data 
DataRn=round(rand(1,NumBits));
DataIn=sign(DataRn);
EbN0=0:1:25;
BER_Interleave=zeros(1,length(EbN0));
BER_Uncoded=zeros(1,length(EbN0));
BER_Coded=zeros(1,length(EbN0));

%------Transmitter side--------    
%Gray mapping of bits into symbols
CmplxBits_Un=length(DataIn)/2;
SigI_Un=zeros(1,CmplxBits_Un);
SigQ_Un=SigI_Un;

%Using the formula for mapping to qpsk
QPSKD_Un=(-2)*sign(DataIn)+1;    
SigI_Un=QPSKD_Un(1:2:NumBits-1);
SigQ_Un=QPSKD_Un(2:2:NumBits);
             
%Complex mapped data
CmplxSg_Un=SigI_Un+1i*SigQ_Un; 

%---Error Correcting Coder------------
%Using single error correcting (15,11) code  
n=15;
k=11;
In=eye(k);     % Idendity Matrix
P=[1 1 1 1;0 1 1 1;1 0 1 1;1 1 0 1;1 1 1 0 ;0 0 1 1;0 1 0 1;...
   0 1 1 0;1 0 1 0;1 0 0 1;1 1 0 0];  % Parity Matrix  
Pt=P.';        % Transpose Matrix
G=[P In];
H=[eye(n-k) Pt];
Ht=H.';         
E=[zeros(1,n);diag(ones(1,n))];  % Error Matrix
S=mod(E*Ht,2);     % Syndrome Table
Se=[S E];          
C_Word=zeros(1,length(DataIn)*n/k);

% Conversion 11 bit messages and code word
C_Word=reshape(DataIn,11,NumBits/11).'; 
Coded_Dt=reshape(mod(C_Word*G,2).',1,(NumBits*15/11)); 

%---------Interleaver-------------------
Depth_In=220; %Matching Depth
CodedDt_In=[];
rg=1;
for rk=1:length(Coded_Dt)/(Depth_In*n)
    Matrix_In=[];
for rep=1:Depth_In
    Matrix_In=vertcat(Matrix_In,Coded_Dt(rg:rg+n-1)  );
    rg=rg+n;
end
CodedDt_In=[CodedDt_In Matrix_In(1:Depth_In*n)];
end

%Gray mapping of coded bits into symbols
CmplxBits_Co=length(Coded_Dt)/2;
SigI_Co=zeros(1,CmplxBits_Co);
SigQ_Co=SigI_Co;

%Using the formula for mapping to qpsk
Coded_Dt=(-2)*sign(Coded_Dt)+1;
SigI_Co=Coded_Dt(1:2:length(Coded_Dt)-1);
SigQ_Co=Coded_Dt(2:2:length(Coded_Dt));

% Complex mapped coded data
CmplxSg_Co=SigI_Co+1i.*SigQ_Co;
    
%Gray mapping of Interleaving bits into symbols
CmplxBits_In=length(CodedDt_In)/2;
SigI_In=zeros(1,CmplxBits_In);
SigQ_In=SigI_In;

%using the formula for mapping to qpsk
CodedDt_In=(-2)*sign(CodedDt_In)+1;
SigI_In=CodedDt_In(1:2:length(CodedDt_In)-1);
SigQ_In=CodedDt_In(2:2:length(CodedDt_In));

% complex mapped data
CmplxSg_In=SigI_In+1i.*SigQ_In;

C=3e8; %velocity of light
fc=10e9; %carrier frequency
lembda=C/fc;
V=(60*1000)/3600; %receiver velocity
fm=V/lembda;
Tc=9/(16*pi*fm); %Coherence Time
Ts=1/NumBits;    %Symbol Time
B2=Tc/Ts;  %block
B=B2.*1.36; %15/11=1.36
B2=round(B2);
B=round(B); 

%noise of fading channel without code
Cnl_Un=zeros(1,length(CmplxSg_Un)+B2);

for nn=1:B2:length(CmplxSg_Un)
        
        Cnl_Un(nn:nn+B2-1)=0.5*(randn*ones(1,B2)+1i*randn*ones(1,B2));
    end
CnlMax_Un=Cnl_Un(1:length(CmplxSg_Un));   

LL=1;
while LL <=length(EbN0)
    %--------CHANNEL-----------
    %Random noise generation and addition to the signal
    SNR=(10)^((EbN0(LL))/10);
    Var=1/sqrt(2*SNR);
    Noise=Var*(randn(1,length(CmplxSg_Un))+ 1i*randn(1,length(CmplxSg_Un))); 
    DataNo_Un=CmplxSg_Un.*CnlMax_Un+Noise;
    
    %--------RECEIVER----------
    %Demapped data and error calculations
    DataNo_Un=DataNo_Un./CnlMax_Un;
    RXBits_Un=zeros(1,NumBits);
   
    %for x=-2*I+1; where the inverse formula is I=-(x-1)/2
    DD=1;
    for nn=1:NumBits/2
       RXBits_Un(DD)= -(sign(real( DataNo_Un(nn)))-1)/2;
       RXBits_Un(DD+1)= -(sign(imag( DataNo_Un(nn)))-1)/2;
       DD=DD+2;
    end
    
     NbEr_In=length(find(DataIn-RXBits_Un));
     BER_Uncoded(LL)=NbEr_In/NumBits;
     LL=LL+1;
end 

semilogy(EbN0,BER_Uncoded,'bo-','LineWidth',1.5,'MarkerEdgeColor','k',...
                 'MarkerSize',5);   
hold on

Cnl_Co=zeros(1,length(CmplxSg_Co)+B);   
    for nnn=1:B:length(CmplxSg_Co)      
        Cnl_Co(nnn:nnn+B-1)=0.5*randn*ones(1,B)+j*0.5*randn*ones(1,B);
    end
    CnlMax_Co=Cnl_Co(1:length(CmplxSg_Co)); 

LLL=1;    
while LLL <= length(EbN0)                   
    %--------CHANNEL-----------
    %Random noise generation and addition to the signal         
    SNR=10.^(EbN0(LLL)/10);
    Var=sqrt((11/15)/(2.*SNR));
    Noise=(Var*randn(1,length(CmplxSg_Co))+ 1i*Var*randn(1,length(CmplxSg_Co))); 
    DataNo_Co=CmplxSg_Co.*CnlMax_Co+Noise ;
    
    %--------RECEIVER----------
    %Demapped data and error calculations
    DataNo_Co=DataNo_Co./CnlMax_Co;
    RXBits_Co=zeros(1,length(Coded_Dt));
                                 % To Decode
    xx=1;
    for jj=1:length(DataNo_Co)
        RXBits_Co(xx)= -(sign(real(DataNo_Co(jj)))-1)/2;
        RXBits_Co(xx+1)= -(sign(imag(DataNo_Co(jj)))-1)/2;
        xx=xx+2;
    end
    %Estimated Code Word
    RXCoded_Co=reshape(RXBits_Co,15,length(Coded_Dt)/15).';
    
    % Syndrome Decoding
    s=mod(RXCoded_Co*Ht,2);
    
    % Error Detection and Correction
    for Row=1:length(Coded_Dt)/15
        for ConT=1:16
            if s(Row,:)==S(ConT,:)
                CorrBits3(Row,:)=xor(RXCoded_Co(Row,:),E(ConT,:));
            end
        end
    end
    
    % Conversion of Message Vector from Code Word
    CorrMsgs3=CorrBits3(:,5:15);
    CorrVect3=reshape(CorrMsgs3',1,NumBits);
    
    % Simulated Bit Error Rate Calculation
    Errors3=find(xor(CorrVect3,DataIn));    
    Errors3=size(Errors3,2);
    BER_Coded(LLL)=Errors3/NumBits;
    LLL=LLL+1; 
end 

semilogy(EbN0,BER_Coded,'gx-','LineWidth',1.5,'MarkerEdgeColor','r',...
                  'MarkerSize',5);  

hold on

Cnl_In=zeros(1,length(CmplxSg_In)+B);   
    for nnn=1:B:length(CmplxSg_In)      
        Cnl_In(nnn:nnn+B-1)=0.5*randn*ones(1,B)+j*0.5*randn*ones(1,B);
    end
    CnlMax_In=Cnl_In(1:length(CmplxSg_In)); 

L=1;
while L <= length(EbN0)                   
    %--------CHANNEL-----------
    %Random noise generation and addition to the signal         
    SNR=10.^(EbN0(L)/10);
    Var=sqrt((11/15)/(2.*SNR));
    Noise=(Var*randn(1,length(CmplxSg_In))+ 1i*Var*randn(1,length(CmplxSg_In))); 
    DataNo_In=CmplxSg_In.*CnlMax_In+Noise ;
    
    %--------RECEIVER----------
    %Demapped data and error calculations
    DataNo_In=DataNo_In./CnlMax_In;
    RXBits_In=zeros(1,length(CodedDt_In));
   
    x=1;
    for j=1:length(DataNo_In)
        RXBits_In(x)= -(sign(real(DataNo_In(j)))-1)/2;
        RXBits_In(x+1)= -(sign(imag(DataNo_In(j)))-1)/2;
        x=x+2;
    end
    
    %---------Deinterleaver-------------------
    Matrix_DeIn=[];
  
    RXCoded_In=[];
    rk=1;
    for rm=1:length(RXBits_In)/(Depth_In*n)
        Matrix_DeIn=[];
    for ri=1:n
        Matrix_DeIn=vertcat(Matrix_DeIn,(RXBits_In(rk:rk+Depth_In-1)));
        rk=rk+Depth_In;    
    end
    RXCoded_In=[RXCoded_In Matrix_DeIn(1:Depth_In*n)];
    end
    MC=1;
    % Syndrome Decoding
    for y=1:n:length(RXCoded_In)    
       S= mod(RXCoded_In(y:y+(n-1))*Ht,2);
       for z=1:n+1
          if S==[Se(z,1) Se(z,2) Se(z,3) Se(z,4)]
             Er=[ Se(z,5) Se(z,6) Se(z,7) Se(z,8) Se(z,9) Se(z,10) Se(z,11) Se(z,12) Se(z,13) Se(z,14) Se(z,15) Se(z,16) Se(z,17) Se(z,18) Se(z,19)];
          end
       end
       CorrCo_In(y:y+(n-1))=xor(RXCoded_In(y:y+(n-1)),Er);
       CorrUn_In(MC:MC+(k-1))=CorrCo_In(y+4:y+(n-1));
       MC=MC+k;
    end    
    % Simulated Bit Error Rate Calculation
    Errors=find(xor(CorrUn_In,DataIn));    
    Errors=size(Errors,2);
    BER_Interleave(L)=Errors/NumBits;
    L=L+1; 
end   

semilogy(EbN0,BER_Interleave,'m^-','LineWidth',1.5,'MarkerEdgeColor','g',...
                  'MarkerSize',5);   
grid on
axis([1 25 10^-5 1]);
xlabel('Eb/No,(dB)')
ylabel('Bit Error Rate')
title('BER Performance Using Interleaver')
legend('UnCoded','coded','Interleaver')