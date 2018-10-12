clc
%bits=111000;
data=randint(1,111000)>0.5; %Generate eqiprobable bits

ebno=0:10;
%BER=zeros(1,length(ebno));
%Formula x=(-2)*data+1 for mapping to qpsk based on consellation diagram
I=(-2)*data(1:2:end)+1;
Q=(-2)*data(2:2:end)+1;
cmplx_data=I+1i*Q;

% Ts=2/Rb, where Rb=1Mb/s and Tc=9/(16*pi*fm),where fm=v/lambda, where lambda=c/fc
% Thus Tc/Ts=160(No.of symbols)

nt1=zeros(1,length(cmplx_data)+160);%Fading Channel Noise

for n=1:160:length(cmplx_data)
    i_ray=0.5*randn*ones(1,160);
    q_ray=0.5*randn*ones(1,160);
    nt1(n:n+159)=(i_ray)+1i*(q_ray);
end
nt=nt1(1:length(cmplx_data));


for i=1:length(ebno) 
    
    % Generation & Adding AWGN to signal
    
    Nd=10.^(ebno(i)/10);
    I_no=1/sqrt(2.*Nd)*randn(1,length(cmplx_data));
    Q_no=1/sqrt(2.*Nd)*randn(1,length(cmplx_data));
    g_noise=I_no+1i*Q_no; 
    
    noisy_sg=cmplx_data.*nt+g_noise;
    
    %Receiver side operation
    rx_data=noisy_sg./nt;
    rx_bits=zeros(1,111000);
    
    % For demapping of bits we will use inverse formula for x=-2*I+1 
    % hence the inverse formula is I=(-x+1)/2
    
    ct=1;
    for ii=1:111000/2
        rx_bits(ct)= -(sign( real( rx_data(ii) ) )-1)/2;
        rx_bits(ct+1)=-(sign( imag( rx_data(ii) ) )-1)/2;
        ct=ct+2;
    end
    
    % Calculation of bit error rate in AWGN channel
    
    err=find(xor(rx_bits,data));    
    err=size(err,2);
    BER(i)=err/111000;
    
end

%For theoretical & practical curve plotting

EbNo=10.^(ebno/10);
theory_BER=0.5*(1-sqrt(EbNo./(1+EbNo)));

pra_BER=.5*erfc(sqrt(10.^(ebno/10)));

% related to graph
semilogy(ebno,BER,'b*-');
hold on
semilogy(ebno,theory_BER,'rs--');
hold on
semilogy(ebno,pra_BER,'g.-');

xlabel('Eb/No(dB)')
ylabel('Bit Error rate(BER)')
title('BER for QPSK-AWGN with Rayleigh fading')
legend('Simulation','Theory','Without Fading')
grid on
    
  



