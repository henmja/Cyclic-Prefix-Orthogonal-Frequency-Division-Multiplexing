clc;clear;close all;

nCh = 10; %Number of multipath channels
nP = 5; %Number of paths per channel
modRate = 10000; %Modulation Rate                
period = 1/modRate;                  
samplingRate = 4;                     
sPeriod = period/samplingRate;                     
cp = 32;

i = 1;
global uwac_params;
uwac_params = cell(nCh,5);

while (i <= nCh) 
    [impRes, prop_delay,processing_gain,ctr,dtr] = UWAC(sPeriod,'path_count',nP);
    if ((isminphase(impRes(:,1)) ~= 0 && size(impRes(:,1),1) < (cp+1)*samplingRate))
        uwac_params{i,1} = impRes;
        uwac_params{i,2} = processing_gain;
        uwac_params{i,3} = prop_delay;
        uwac_params{i,4} = ctr;
        uwac_params{i,5} = dtr;
        disp(i);
        i = i+1;
    end
end

MCSimulations = 5;
UWACplot(MCSimulations, uwac_params);

function [] = UWACplot(numSims, uwac_params)
MCsimulations = numSims;
Hn = uwac_params{MCsimulations,1}; %Finite impulse response filter (tap vector = H0,H1,..,Hn)
processingGain  = uwac_params{MCsimulations,2};
propagationDelay = uwac_params{MCsimulations,3}; 

figure(1);
name = 'Sans Serif';
size = 15;
plot(Hn);
xlabel('$Time (Sample)$','fontname',name,'fontsize',size,'interpreter','latex');
ylabel('$Impulse Response$','fontname',name,'fontsize',size,'interpreter','latex');

figure(2);
%processing gain y and propagation delay x (the gains in the plot are the average gains for all 
%the paths):
bar(propagationDelay(:,1)*1000,processingGain);
xlabel('$Delay$','fontname',name,'fontsize',size,'interpreter','latex');
ylabel('$Gain$','fontname',name,'fontsize',size,'interpreter','latex');
set(gca,'fontname',name,'fontsize',size);
end

function [impulseResponse,lag,taps,varargout] = UWAC(period,varargin)
p = inputParser;
addParameter(p,'guardInterval',0.0246);
addParameter(p,'path_count',15);
addParameter(p,'lag_mean',0.001);
addParameter(p,'powerChange',20);
addParameter(p,'velocity',20);
addRequired(p,'period',@isnumeric);
parse(p, period, varargin{:});

period = p.Results.period;
path_count = p.Results.path_count;
lag_mean = p.Results.lag_mean;
powerChange = p.Results.powerChange;
guardInterval = p.Results.guardInterval;
velocity = p.Results.velocity;

%Doppler frequency:
c = 1500; %velocity of the waves in the medium
f = velocity/c;
f_n = length(f);
[ctr,dtr] = rat(f+1);
varargout{1} = ctr;
varargout{2} = dtr;

%Rayleigh Multipath:
pd = makedist('exponential','mu',lag_mean);
R     = random(pd,path_count,1);
Y = ceil(R/period);
lag_location = cumsum(Y);
lag = period*(lag_location-1);
a = log(10^(powerChange/10))/guardInterval;
avgPower = exp(lag*-a);
taps = raylrnd(avgPower*sqrt(2/pi));

B = repmat(f',path_count,1);
temp = repmat(lag_location,1,f_n);
new_lag_loc = ceil(temp/(1+B));
impulseResponse = zeros(max(max(new_lag_loc)),f_n);

for i = 1:f_n
    impulseResponse(new_lag_loc(:,i),i) = taps;
    impulseResponse(:,i) = impulseResponse(:,i)/norm(impulseResponse(:,i));
end
lag = (new_lag_loc - 1)*period;
end
