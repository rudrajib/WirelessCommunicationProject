clc
clear all
bits=100000;
data=randint(1,bits)>0.5;
%---debugging---
%data=[1 1 1]
%xxxxxxxxxx
ebno=0:15;
BER=zeros(1,length(ebno));
   
    %---Transmitter---------
    %Gray mapping of bits into symbols
    col=length(data)/2;
    I=zeros(1,col);
    Q=I;
    
    I=data(1:2:bits-1);
    Q=data(2:2:bits);
    
    I= -2.*I+1;
    Q= -2.*Q+1;
    
    symb=I+j.*Q;
    
    %noise of fading channel
    htt1=zeros(1,length(symb)+160);
    
   
    for n=1:160:length(symb)
        
        htt1(n:n+159)=0.5*randn*ones(1,160)+1i*0.5*randn*ones(1,160);
    end
    ht=htt1(1:length(symb));        
    
    
          
for i=1:length(ebno)
           
    
    %--------CHANNEL-----------
    %Random noise generation and addition to the signal
    npsd=10.^(ebno(i)/10);
    n_var=1/sqrt(2.*npsd);
    gaussian_noise=n_var*randn(1,length(symb))  +1i*n_var*randn(1,length(symb));
    
    
    rx_symb=symb.*ht+gaussian_noise;
    %xxxxxxxxxxxxxxxxxxxxxxxxxx
    
    
    
    %-------RECEIVER-----------
    
    rx_symb=rx_symb./ht;
    
    recv_bits=zeros(1,bits);
    %demapping
    k=1;
    for ii=1:bits/2
        recv_bits(k)=  -( sign(  real(  rx_symb(ii)  )  )  -1)/2;
        recv_bits(k+1)=-( sign(  imag(  rx_symb(ii)  )  )  -1)/2;
        k=k+2;
    end
        
           
   %---SIMULATED BIT ERROR RATE----
    errors=find(xor(recv_bits,data));    
    errors=size(errors,2);
    BER(i)=errors/bits;
    %xxxxxxxxxxxxxxxxxxxxxxxxxxx
end

semilogy(ebno,BER,'b.-');
hold on

ebn=10.^(ebno/10);
thr=.5.*(1-sqrt(ebn./(1+ebn)));

pb=.5*erfc(sqrt(10.^(ebno/10)));
semilogy(ebno,thr,'rx-');
hold on
semilogy(ebno,pb,'g+-');
xlabel('Eb/No (dB)')
ylabel('Bit Error rate')
title('Performance of QPSK with Rayleigh Fading')
legend('Simulation','Theory','Without Fading')
grid on