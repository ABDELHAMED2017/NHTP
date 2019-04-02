% demon sparse logistic regression problems with real data
clc; close all; clear;

prob      = 'newsgroup'; %'colon-cancer'
measure   = load(strcat(prob,'.mat')); 
label     = load(strcat(prob,'_label.mat'));   
label.b(label.b==-1)= 0;
[m,n]     = size(measure.A);

normtype  = 1;    
if m     >= 1000; normtype=2; end

data.A    = normalization(measure.A,normtype); 
data.At   = data.A';
data.b    = label.b; 
s         = ceil(0.05*m);
pars.Draw = 1; 
out       = NHTP(n,s,data,'SLR',pars) 

if isfield(pars,'Draw') && pars.Draw
saveas(figure(1), [pwd strcat(strcat('/outputs/',prob))]);  
saveas(figure(1), [pwd strcat(strcat('/outputs/',prob),'.eps')],'epsc');
end
