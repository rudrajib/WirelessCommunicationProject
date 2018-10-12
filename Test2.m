%--------------------------------------------------------------------------
%In the second part of the project we have to simulate the communication
% system having Doppler Spread due to mobile movement in the channel.we 
%Simulate and draw the bit error rate of the system for maximum mobile
%speed of 81 Km/Hour. We Assumed that the channel is flat and no equalizer
%is used in the system
%--------------------------------------------------------------------------
clc
clear all;
close all;
datalength=2e5;             % Total number of  data for transmission
Eb_N0_dB=0:3:45;                 % multiple Eb/N0 values


%--------------------------------------------------------------------------
    % As Tc=9/(16*pi*fm),where fm=v/lembda, where lembda=c/fc, where
    % Ts=2/Rb(for the case of QPSK), where Rb=1Mb/s, so Tc/Ts=no. of symbol
    % so we have total no of symbols are 489.4 per block so we round to 500
%--------------------------------------------------------------------------
    symbolsperblock=500;           % no.of symbol per block
    total_sy=datalength/2;
    totalblocks=total_sy/symbolsperblock; % tolal no. of blocks
     data=round(rand(1,datalength)); % generating 0,1 with equal probability
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %using the formula of x=-2*data+1 for mapping to qpsk
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    i_data = data(1:2:end);
    q_data = data(2:2:end);
    I_qpsk = (-2)*i_data+1;
    Q_qpsk = (-2)*q_data+1;
    
      symbol1=I_qpsk+sqrt(-1)*Q_qpsk;
      
   %-----------------------------------------------------------------------
   % Co-efficient of rayleigh fading%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   I_ray=sqrt(0.5)*(randn(1,totalblocks));
   Q_ray=sqrt(0.5)*(randn(1,totalblocks));
   coefficient_ray=I_ray+sqrt(-1)*Q_ray;
   %multiply rayleigh co-efficient to signal
%--------------------------------------------------------------------------
for rs = 1:length(Eb_N0_dB)
  
   
   %converting to 1000 by 500 matrix
   symbolMatrix = reshape(symbol1,totalblocks,symbolsperblock);
   
   counter=1;
   for t=1:totalblocks
           
   RalyMatrix(t,:)= symbolMatrix(t,:)*coefficient_ray(t);
   counter=counter+1;
   end
   %-----------------------------------------------------------------------
   %%%%%%%%%%%%%%% generation of  AWGN %%%%%%%%%%%%%%%%%
  
   
   No=(10.^(-Eb_N0_dB(rs)/10));
   I_noise=randn(1,datalength/2)*(sqrt(No/2));
   Q_noise=randn(1,datalength/2)*(sqrt(No/2));
   noise=I_noise+sqrt(-1)*Q_noise;
   
   noiseMatrix = reshape(noise,totalblocks,symbolsperblock);
  
   % Add noise in signal
   transmitted_data = RalyMatrix + noiseMatrix;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % At the receiver side, data demapping and error calculations
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   counter1=1;
   for k=1:totalblocks
           received_data(k,:)=transmitted_data(k,:)/coefficient_ray(k);
           counter1=counter1+1;

   end
   
   received_data1 = reshape(received_data, 1,datalength/2);
   mat=[real(received_data1); imag(received_data1)];
   vector=reshape(mat,1,datalength);
   %-----------------------------------------------------------------------
   % get demapped bits from the receive signal using inverse formula 
   %for x=-2*I+1; where the inverse formula is I=(-x+1)/2%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   vector1=sign(vector);
   received_bits1=(-vector1+1)/2;
   received_bits=round(received_bits1);
   ray_error(rs)=sum (xor(data,received_bits))/datalength;
    
    %----------------------------------------------------------------------
    % for the purpose of comparision to that of the rayleigh model we
    % also want to find the BER IN AWGN CHANNEL and then compare both
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    In_phase=(-2)*sign(data(1:2:end))+1;
    Q_phase=(-2)*sign(data(2:2:end))+1;
    symbol2=In_phase+sqrt(-1)*Q_phase;
    % Channel - AWGN
    N1o=(10.^(-Eb_N0_dB(rs)/10));
    awgn_I=randn(1,length(data)/2)*sqrt(N1o/2);
    awgn_Q=randn(1,length(data)/2)*sqrt(N1o/2);
    noise1=awgn_I+sqrt(-1)*awgn_Q;
    trans_data=symbol2+noise1;
    pk=[real(trans_data) ;imag(trans_data) ];
    x2=reshape(pk,1,datalength);
    est_data1=round((sign(x2)-1)/(-2));
    
 %  calculate errror without rayleigh
   errors =sum(xor(data,est_data1));
   Ber_awgn(rs)=errors/datalength;
   
   
end
 %-------------------------------------------------------------------------
 % genaration of the theoritical curve for Rayleigh fading 
 

EbN0 = 10.^(Eb_N0_dB/10);
theoryBerRay = 0.5*(1-sqrt(EbN0./(EbN0+1)));
figure
semilogy(Eb_N0_dB,theoryBerRay,'rs-','LineWidth',4,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','y',...
                'MarkerSize',4);
hold on
semilogy(Eb_N0_dB,ray_error,'g*-','LineWidth',2,...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','r',...
                'MarkerSize',2);
            hold on;
            semilogy(Eb_N0_dB,Ber_awgn,'r*-','LineWidth',3,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',8);
            grid on;
legend('Theoritically BER in Rayleigh Fading channel',...
         'Simulation,BER in Rayleigh fading channel','BER in AWGN channel');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('BER for QPSK in AWGN Channel and Rayleigh fading Channel');