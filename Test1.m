%--------------------------------------------------------------------------
%This Program is the simulation of the communication system to draw the 
%charachteristic curve of the bit error rate of the simulated data versus
%Eb/No(dB) and compare with the bit error rate of the theoritical data
%(Characteristic curve of qpsk)versus Eb/No(dB)
clc;
clear all;
close all;
datalength=1e6; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generating random data and round to one's and zeros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datarandom=round(rand(1,datalength));
data=sign(datarandom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%using the formula of x=-2*data+1 for mapping to qpsk
%according to constellation diagram given below

qpsk_data=(-2)*sign(data)+1;              
I=qpsk_data(1:2:end);
Q=qpsk_data(2:2:end);
% complex mapped data
complex_data=I+1i*Q;                 
%         ^
%         |
%  10 x   |   x 00   (odd bit, even bit)
%         |
%  -------+------->  real part (I channel)
%         |
%  11 x   |   x 01
%         |
% constellation diagram of the mapped data
scatterplot(complex_data); 
title(' constellation diagram for transmitted data having no Noise and qpsk mapped ');
hold on ;
grid on;
EbNodB=0:2:10;
for i=1:length(EbNodB);
%generation of the AWGN 
   No=(10)^((-EbNodB(i))/10);
    AWGN= randn(1,datalength)*sqrt(No/2);
    Ni= AWGN(1:2:end);
    Nq=AWGN(2:2:end) ;
    %Complex noise
    AWGN_Complex=Ni+1i*Nq; 
    
    Noisy_signal=complex_data+ AWGN_Complex;
    
    %----------------------------------------------------------------------
    % At the receiver side demapped data and error calculations
    
    g=[real(Noisy_signal);imag(Noisy_signal)];
     h=reshape(g,1,datalength); 
     k= sign(h);
   
    % get demapped bits from the receive signal using inverse formula 
    %for x=-2*I+1; where the inverse formula is I=(-x+1)/2

     m=(-k+1)/2;
    demapped_data=round(m);
    no_of_biterror=length(find(xor(datarandom,demapped_data)));
    BER(i)=no_of_biterror/datalength ; 
end
%--------------------------------------------------------------------------
% constellation diaram of the noisy signal
scatterplot(Noisy_signal);
title('Constellation diagram for signal after adding noise');
hold on;
grid on;
% formula for the calculation of theoritical Bit error rate
theoretical_ber = 0.5*erfc(sqrt(10.^(EbNodB/10)));

%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%This is all about graph
figure
semilogy(EbNodB,theoretical_ber,'rs--','LineWidth',4,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',10);
hold on
semilogy(EbNodB,BER,'g*-','LineWidth',2,...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','r',...
                'MarkerSize',10);
grid on;
legend('theoretical', 'simulation');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER for QPSK in AWGN channel');