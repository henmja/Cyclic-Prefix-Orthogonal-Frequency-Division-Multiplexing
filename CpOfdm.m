clear;close all;clc;


% ------ Remember to run UWACmain.m before running this file, in order to
% generate the global UWAC parameters.


global uwac_params; 
disp(uwac_params);

SNR  = 0:2:30; %Signal to noise ratio in decibels
SNRlen = length(SNR); %size of SNR

bitsPerSymbol = 4; %16QAM = 2^4QAM -> 4 bits per symbol.
constellationSize = 2^bitsPerSymbol;


nGB = [16; 16]; %number of guard band subcarriers allocated to the left and right guard bands.
nFFT = 1024; %total number of subcarriers.
nBlock = nFFT-(sum(nGB)); %Subcarriers per block.


%Parameters for encoding:
K = 3; %Constraint Length
tableLength = 5*K; 
x = 5; 
k = 7;
codeRate = 0.5; %Convolutional Code Rate.
trellis = poly2trellis(K,[x k]);

generatedBits = bitsPerSymbol*(nBlock*codeRate);

interleaver = randi([1 1e6],1,1);

symbolRate = 24*1000;
symbolPeriod = 1/symbolRate;
BW = 5*1000;

%Convert sampled digital signal to a that of a higher sampling rate.
UpSampling = ceil(symbolRate/BW);
samplingPeriod = symbolPeriod/UpSampling;
samplingFrequency = 1/samplingPeriod;

carrierFrequency = 15000;


% Low-Pass-Filter (LPF)
filterOrder = 1000; %Maximum amount of delay to calculate output
cutoffFrequency = 0.5*samplingFrequency/UpSampling; %frequencies below pass, others attenuate.
angularFrequency = 2*pi*cutoffFrequency/samplingFrequency; %rotation rate
x = 1:filterOrder/2;
tempFilter = angularFrequency/pi; %convolution
tempFilter2 = sin(angularFrequency*x)./(pi*x); %convolution's impulse response
LPF = [fliplr(tempFilter2) tempFilter tempFilter2]; %Low-Pass-Filter
LPF_Hamming = LPF.*hamming(filterOrder+1).'; %Hamming window (0 outside chosen interval).

if((mod(filterOrder,2) == 0))
    transientDelay = (filterOrder)/2; %time it takes response to reach half the final value.
else
    transientDelay = (1+filterOrder)/2;
end

nMC =  1; %Monte carlo iterations
BER = zeros(nMC,SNRlen); %Bit Error Rate init





for i = 1:nMC
 
    y = uwac_params{i,1}; %Set Channel
    %minimum equal to length of multipath channel:
    CyclicPrefixLength = ceil(size(y,1)/UpSampling) - 1;
    
    bits = randi([0 1],generatedBits,1);
    
    RCPMat = [zeros(nFFT,CyclicPrefixLength), eye(nFFT)]; %Receiver's matrix.
    
    TCode = convenc(bits,trellis); %transmitters code.
    
    Tintrlv = randintrlv(TCode,interleaver); %interleave/randomly rearrange TCode.
    
    %Modulate interleaved code:
    Tmod = qammod(Tintrlv,constellationSize,'gray','InputType','bit','UnitAveragePower',true);
       

    symbols = [zeros(nGB(1),1); Tmod; zeros(nGB(2),1)];
    T_IFFT = sqrt(nBlock)*ifft(ifftshift(symbols)); %OFDM Modulation - IFFT
    
    T_CP = [T_IFFT(end-CyclicPrefixLength+1:end,:); T_IFFT].'; %prepend CP
    
    T_CP_Up = upsample(T_CP,UpSampling); %Convert from digital to analog
    
    T_Filter = fftfilt(LPF_Hamming,T_CP_Up);    
    
    
    timeVector = samplingPeriod*(0:(size(T_Filter,2)-1));
    carrier = exp(2*1i*pi*carrierFrequency*timeVector); %carrier at frequency, carrierfrequency.
    T_Pass = real(carrier.*T_Filter); % Convert baseband signal to passband signal
    T_Pwr = norm(T_Pass)^2/length(T_Pass); %Mean transmission power
    
    R_Ch = fftfilt(y(:,1),T_Pass); %receiver channel
    
    %receiver side
    
    for j = 1:SNRlen
        
        %AWGN (Add White Gaussian Noise)
        noiseVector  = randn(1,size(R_Ch,2)); %noise
        % Average noise power
        avgPower = norm(noiseVector)^2/length(noiseVector);
        R_Ch_AWGN = R_Ch + sqrt(((T_Pwr/avgPower))*(10^(-SNR(j)/10)))*noiseVector;
        
        R_Baseband = (R_Ch_AWGN.*conj(carrier)).'; %Demodulate
        R_Temp = fftfilt(LPF_Hamming,[R_Baseband; zeros(filterOrder-1,1)]);
        R_Temp2 = R_Temp(transientDelay+1:end-transientDelay);
        R_DownSample = downsample(R_Temp2,UpSampling);%Reverse Pulse Shaping
        
        %Estimating channel information for OFDM:
        pilotVector = conj(symbols)./(abs(symbols).^2+10^(-SNR(j)/10));
        estVector = fftshift(fft(RCPMat*R_DownSample,nFFT)).*pilotVector/sqrt(nFFT);
        %Equalization matrix:
        estMat = diag(conj(estVector)./(abs(estVector).^2+10^(-SNR(j)/10)));
        
        R_EstSymbTemp = estMat*fftshift(fft(RCPMat*R_DownSample,nFFT))/sqrt(nFFT);
        R_EstSymb = R_EstSymbTemp(nGB(1)+1:end-nGB(2)); %Remove null subcarries

        %Digital demodulation:
        R_Demodulated = qamdemod(R_EstSymb,constellationSize,'gray','OutputType',...
            'approxllr','UnitAveragePower',true,'NoiseVariance',avgPower);
  
        % Deinterleaving
        R_Deintrlv = randdeintrlv(R_Demodulated,interleaver);
 
        %Convolution decoding using the Viterbi algorithm 
        estBits = vitdec(R_Deintrlv,trellis,tableLength,'cont','unquant');
        
        %Estimated Bit Error rate:
        BER(i,j) = biterr(bits(1:end-tableLength),estBits(tableLength+1:end))/(generatedBits);
    end   
end

avgBER = mean(BER,1);
disp(avgBER);
disp(SNR);
figure(1);
%Signal to noise ratio and BER
plot(SNR,avgBER);
name = 'Sans Serif';
size = 15;
xlabel('$SNR$','fontname',name,'fontsize',size,'interpreter','latex');
ylabel('$BER$','fontname',name,'fontsize',size,'interpreter','latex');