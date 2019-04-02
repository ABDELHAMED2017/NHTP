% demon compressed sensing problems with randomly generated data
clc; clear; close all;

n       = 2000;  
m       = ceil(n/4); 
s       = ceil(0.05*n);                      
 
% % You could input any data including (data.A, data.At, data.b) 
I         = randperm(n); 
x         = zeros(n,1);  
x(I(1:s)) = randn(s,1);
data.x_opt= x;
data.A    = randn(m,n)/sqrt(m);
data.At   = data.A';
data.b    = data.A*data.x_opt;  

% Or you could input data from our data generation function
% ExMat = 1;
% MatType = {'GaussianMat','PartialDCTMat'}; 
% data    = compressed_sensing_data(MatType{ExMat}, m,n,s,0);

pars.Draw = 1;
out       = NHTP(n,s,data,'CS',pars); 

fprintf('\n Sample size:       m=%d,n=%d\n', m,n);
fprintf(' Recovery time:    %6.3fsec\n',  out.time);
if isfield(data,'x_opt')
fprintf(' Recovery accuracy: %5.2e\n\n', norm(out.sol-data.x_opt)/norm(data.x_opt));
end