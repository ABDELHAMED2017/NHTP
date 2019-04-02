clc; clear; close all;

ExNam   = 1; %= 1, 2, 3, 4 or 5                          
ExMat   = 1; %= 1 or 2

n       = 1000;  
m       = ceil(0.25*n); 
s       = ceil(0.05*n);
if ExNam> 3
n       = 2; 
s       = 1;
end   
    
switch ExNam
    
    case 1 % demon compressed sensing problems
    MatType = {'GaussianMat','PartialDCTMat'}; 
    data    = compressed_sensing_data(MatType{ExMat}, m,n,s,0);
    probname='CS'; 
    
    case 2 % demon sparse logistic regression problems
    MatType = {'Indipendent','Correlated'};
    rho     = 0.5; % 0<= rho <=1 it is useful for 'Correlated' data
    data    = logistic_random_data(MatType{ExMat},m,n,s,rho);
    probname='SLR';   
    
    case 3 % demon a simple example 
    a       = 0.1*randn; b=0.1*rand(n,1);
    data    = @(var,flag)simple_ex2_func(var,flag, a, b);
    probname='SCO'; 
    
    case 4 % demon a simple example 
    data    = @(var,flag)simple_ex3_func(var,flag);
    probname='SCO';  
    
    case 5 % demon a simple example 
    data    = @(var,flag)simple_ex4_func(var,flag);
    probname='SCO';  
end
 
pars.Draw = 1;
out       = NHTP(n,s,data,probname,pars)
 