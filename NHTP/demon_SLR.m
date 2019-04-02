% demon sparse logistic regression problems with randomly generated data
clc; clear; close all;

n    = 10000;  
m    = ceil(n/5); 
s    = ceil(0.05*n);
                      
% % You could input any data including (data.A, data.At, data.b) 
% I0      = randperm(n); 
% I       = I0(1:s);
% x       = zeros(n,1);
% x(I)    = randn(s,1);
% data.A  = randn(m,n); 
% data.At = data.A';
% q       = 1./(1+exp(-data.A(:,I)*x(I)));
% data.b  = zeros(m,1);
% for i   = 1:m; data.b(i) = randsrc(1,1,[0 1; 1-q(i) q(i)]); end 

% Or you could input data from our data generation function 
ExMat = 1;
MatType = {'Indipendent','Correlated'};
data    = logistic_random_data(MatType{ExMat},m,n,s,0.5);
   
pars.Draw = 1;
out       = NHTP(n,s,data,'SLR',pars);

fprintf('\n Sample size:   m=%4d,n=%4d\n', m,n);
fprintf(' CPU time:     %6.3fsec\n',  out.time);
fprintf(' Logistic Loss: %5.2e\n\n', out.obj);