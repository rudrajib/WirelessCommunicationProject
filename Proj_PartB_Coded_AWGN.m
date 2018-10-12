clc;
clear;
NumBits=220000;              % Number of bits 
DataIn=round(rand(1,NumBits));   % To generate input bits of length NumBits, later must divide it to imaginary I and Real Q
Min=0;
Max=10;
EbN0=(Min:Max);  % Range of signal to Noise ratio, Eb/N0, horizontally  
BER=zeros(1,length(EbN0));   % BER, Bit Error Rate
                       %   Transmitter side
              % Coding, single error correcting (15,11)code where n=15, k=11 
n=15;
k=11;
In=eye(k);     % Idendity Matrix In
               % Parity Matrix P 
P=[1 1 1 1;0 1 1 1;1 0 1 1;1 1 0 1;1 1 1 0 ;0 0 1 1;0 1 0 1;0 1 1 0;1 0 1 0;1 0 0 1;1 1 0 0];
Pt=P.';         % Pt is transpose matrix of matrix P
G=[P In];
H=[eye(n-k) Pt];
Ht=H.';         %  Ht is transpose matrix of Matrix H
E=[zeros(1,n);diag(ones(1,n))];  % Error Matrix E
Se=[mod(E*Ht,2) E];          % Syndrome table Se
W=zeros(1,length(DataIn)*n/k);
M=1;
for i=1:k:length(DataIn)
    W(M:M+(n-1))=mod(DataIn(i:i+(k-1))*G,2);
    M=M+n;
end                                  
CodedData=W;      % Coded data                     
                            % Transmitter side                              
                       % Gray mapping of bits into symbols
Co=length(CodedData)/2;
SigI=CodedData(1:2:length(CodedData)-1);
SigQ=CodedData(2:2:length(CodedData));
SigI= -2.*SigI+1;
SigQ= -2.*SigQ+1;
symbol=SigI+1i.*SigQ;
L=1;
while L <= length(EbN0)                   
              % AWGN Addetive White Guassian Noise,to generate Random Noise and add to complex signal
    NoisS=10.^(EbN0(L)/10);
    NoisVariable=1/sqrt(2.*NoisS);
    RxSymbol=symbol+(NoisVariable*randn(1,length(symbol))+ 1i*NoisVariable*randn(1,length(symbol)));           
                              %  At the Reciver side:
    RX_Bits=zeros(1,length(CodedData));
                                 % To Decode
    x=1;
    for j=1:length(RxSymbol)
        RX_Bits(x)= -(sign(real(RxSymbol(j)))-1)/2;
        RX_Bits(x+1)= -(sign(imag(RxSymbol(j)))-1)/2;
        x=x+2;
    end
    RX_CodedBits=RX_Bits;                            
    CodedBitsCorected=zeros(1,length(CodedData));
    BitsUnCoded=zeros(1,NumBits);
    length(DataIn);
    MC=1;
    for y=1:n:length(RX_CodedBits)    
       S= mod(RX_CodedBits(y:y+(n-1))*Ht,2);
       for z=1:n+1
          if S==[Se(z,1) Se(z,2) Se(z,3) Se(z,4)]
             Er=[ Se(z,5) Se(z,6) Se(z,7) Se(z,8) Se(z,9) Se(z,10) Se(z,11) Se(z,12) Se(z,13) Se(z,14) Se(z,15) Se(z,16) Se(z,17) Se(z,18) Se(z,19)];
          end
       end
       CodedBitsCorected(y:y+(n-1))=xor(RX_CodedBits(y:y+(n-1)),Er);
       BitsUnCoded(MC:MC+(k-1))=CodedBitsCorected(y+4:y+(n-1));
       MC=MC+k;
    end
    DataIn;
    W;
                          % CodedData;
    RX_CodedBits;
    RX2_Bits=BitsUnCoded;       
                          %  BER, Bit Error Rate by Simulation 
    Errors=find(xor(RX2_Bits,DataIn));    
    Errors=size(Errors,2);
    BER(L)=Errors/NumBits;
    L=L+1; 
end          % end of while loop
             % Using fotmula to calculate Theroretical Bit Error Rate values Coded
                             
Theory_BER =0.5*erfc(sqrt(10.^(EbN0/10)));  
               % Using the the Sumation formula to find PB, Uncoded
t=1;
Theory_PB=zeros(1,length(EbN0));
i=1;
while i<=length(EbN0)
    for j=t+1:n
        Theory_PB(i)=Theory_PB(i)+j*nchoosek(n,j)*Theory_BER(i).^j*(1-Theory_BER(i)).^(n-j);
    end
    Theory_PB(i)=Theory_PB(i)/n;
    i=i+1;
end        % end of while loop
                        % To plot and lable on the semilog graphs
semilogy(EbN0,BER,'bx-');    % Bit Error Rate, BER Simulation
hold on                 
semilogy(EbN0,Theory_BER,'r^-');  % Theoretical graph, Uncoded
hold on
semilogy(EbN0,Theory_PB,'mo-');    % Theoretical graph, Coded       
xlabel('Eb/No,(dB)')
ylabel('Bit Error Rate')
title('QPSK Performance with coded and AWGN')
legend('Simulation','Theoretical, UnCoded','Theoretical, Coded')
grid on