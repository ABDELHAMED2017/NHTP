% This code presents the success recovery rate of NHTP
clc; clear; close all; 

test       = 1; %=1 succ rate v.s. s; =2 succ rate v.s. m/n
ExMat      = 1; %=1 Gaussian matrix;  =2 Partial DCT matrix

n           = 256; 
m           = floor(0.25*n);
s           = floor(0.05*n);
noS         = 500;
MatType     = {'GaussianMat','PartialDCTMat'};
if test  == 1
test0       = floor(0.16*m):2:floor(0.6*m); 
else
test0       = 0.08:0.02:0.22;
end
 
SuccRate    = [];
pars.IterOn = 0;
pars.Draw   = 0;
 
for j=1:length(test0) 
    rate = 0; 
    if test==1  
        s = test0(j);
    else           
        m = floor(test0(j)*n);
    end    
    for S = 1:noS         
        data  = compressed_sensing_data(MatType{ExMat},m,n,s,0 );       
        out   = NHTP(n,s,data,'CS',pars);  clc; SuccRate     
        if norm(out.sol-data.x_opt)/norm(data.x_opt)<1e-2; rate=rate+1;end 
    end
    clc; SuccRate  = [SuccRate rate]    
end

figure
plot(test0,SuccRate/noS,'r*-') ; hold on
if test==1; xlabel('s')  ;
else;          xlabel('m/n');
end
ylabel('Success Rate');
axis([min(test0) max(test0) 0 1]); grid on;
legend('NHTP','Location','NorthEast'); hold on 

