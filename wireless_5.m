clc
clear all
bits=242000;
gen_data=randint(1,bits)>0.5;
%gen_data=[1 1 1 1 0 0 1 1 0 1 1]
%---debugging---
%data=[1 1 1]
%xxxxxxxxxx
ebno=0:25;
BER=zeros(1,length(ebno));
   
    %---Transmitter---------
    
    %coding
n=15;
k=11;
P=[1 1 1 1 ;0 1 1 1 ;1 0 1 1;1 1 0 1;1 1 1 0;0 0 1 1;...
   0 1 0 1 ;0 1 1 0 ;1 0 1 0;1 0 0 1;...
   1 1 0 0];
G=[P eye(k)];
H=[eye(n-k) P.'];
Ht=H.';
e=[zeros(1,n);diag(ones(1,n))];
synd_table=[ mod(e*Ht,2) e];
 
U=zeros(1,length(gen_data)*n/k);
kk=1;
for ii=1:k:length(gen_data)
    U(kk:kk+(n-1))=mod(gen_data(ii:ii+(k-1))*G,2);
    kk=kk+n;
end
 
%xxxxxxxx

       
    
 
%coded data
U;

%interleaver

depth=220;

U2=[];
rg=1;
for rk=1:length(U)/(depth*n)
    inter_mat=[];
for rep=1:depth
    inter_mat=vertcat(inter_mat,U(rg:rg+n-1)  );
    rg=rg+n;
end

U2=[U2 inter_mat(1:depth*n)];
end
U=[];
U=U2;


data=U;
    
   
    %Gray mapping of bits into symbols
    col=length(data)/2;
    I=zeros(1,col);
    Q=I;
    
    I=data(1:2:length(data)-1);
    Q=data(2:2:length(data));
    
    I= -2.*I+1;
    Q= -2.*Q+1;
    
    symb=I+j.*Q;
    
            %noise of fading channel
    htt1=zeros(1,length(symb)+160);
    
   
    for n1=1:160:length(symb)
        
        htt1(n1:n1+159)=0.5*randn*ones(1,160)+j*0.5*randn*ones(1,160);
    end
    ht=htt1(1:length(symb)); 
          
for i=1:length(ebno)
           
    
    %--------CHANNEL-----------
    %Random noise generation and addition to the signal
    npsd=10.^(ebno(i)/10);
    n_var=sqrt((k/n)./(2.*npsd));
    gaussian_noise=n_var*randn(1,length(symb))+j*n_var*randn(1,length(symb));
    rx_symb=symb.*ht+gaussian_noise ;
   
    %xxxxxxxxxxxxxxxxxxxxxxxxxx
    
    
    %-------RECEIVER-----------
        rx_symb=rx_symb./ht;
        %rx_symb=rx_symb;
    recv_bits=zeros(1,length(data));
    
    %demapping
    kt=1;
    for ii=1:length(rx_symb)
        recv_bits(kt)  =  -( sign(  real(  rx_symb(ii)  )  )  -1)/2;
        recv_bits(kt+1)=  -( sign(  imag(  rx_symb(ii)  )  )  -1)/2;
        kt=kt+2;
    end

    %recvd_coded_bits=recv_bits;
    %decoding

deinter_mat=[];
  
    recvd_coded_bits=[];
    rk=1;
    for rm=1:length(recv_bits)/(depth*n)
        deinter_mat=[];
    for ri=1:n
        deinter_mat=vertcat(deinter_mat,(recv_bits(rk:rk+depth-1)));
        rk=rk+depth;
        
    end
    recvd_coded_bits=[recvd_coded_bits deinter_mat(1:depth*n)];
    end

corrected_coded_bits=zeros(1,length(U));
uncoded_bits=zeros(1,bits);
length(gen_data);
 
c=1;
for kkk=1:n:length(recvd_coded_bits)    
S=mod(recvd_coded_bits(kkk:kkk+(n-1))*Ht,2);
for iii=1:n+1
    if S==[synd_table(iii,1) synd_table(iii,2) synd_table(iii,3) synd_table(iii,4)]
        
        ed=[ synd_table(iii,5) synd_table(iii,6) ...
            synd_table(iii,7) synd_table(iii,8) synd_table(iii,9) ...
            synd_table(iii,10) synd_table(iii,11) synd_table(iii,12) ...
            synd_table(iii,13) synd_table(iii,14) synd_table(iii,15)...
             synd_table(iii,16) synd_table(iii,17) synd_table(iii,18) synd_table(iii,19)];
    end
end
corrected_coded_bits(kkk:kkk+(n-1))=xor(recvd_coded_bits(kkk:kkk+(n-1)),ed);

uncoded_bits(c:c+(k-1))=corrected_coded_bits(kkk+4:kkk+(n-1));

c=c+k;
end
    gen_data;
    U;
    recvd_coded_bits;
    recv_bits2=uncoded_bits;
    
    
           
   %---SIMULATED BIT ERROR RATE----
    errors=find(xor(recv_bits2,gen_data));    
    errors=size(errors,2);
    BER(i)=errors/bits;
    %xxxxxxxxxxxxxxxxxxxxxxxxxxx
end

figure

semilogy(ebno,BER,'b.-');
hold on

ebn=(10.^(ebno/10));
ebn=ebn.*(k/n);
pc=.5.*(1-sqrt(ebn./(1+ebn)));

t=1;
Pb=zeros(1,length(ebno));
for i=1:length(ebn)
    for j=t+1:n
        Pb(i)=Pb(i)+j*nchoosek(n,j)*pc(i).^j*(1-pc(i)).^(n-j);
    end
    Pb(i)=Pb(i)./n;
end


semilogy(ebno,pc,'rx-');

%theoritical performance without fading

pc=0.5*erfc(sqrt(ebn.*(11/15)));

t=1;
Pb=zeros(1,length(ebno));
for i=1:length(ebno)
    for j=t+1:n
        Pb(i)=Pb(i)+j*nchoosek(n,j)*pc(i).^j*(1-pc(i)).^(n-j);
    end
    Pb(i)=Pb(i)./n;
end
hold on
semilogy(ebno(1:10),Pb(1:10),'g+-');

xlabel('Eb/No (dB)')
ylabel('Bit Error rate')
title('Performance of Coded QPSK with Interleaver in Rayleigh Fading')
legend('Simulation','Theory','Without Fading')
grid on
