function [out] = simple_ex3_func(x,fgH)

% This code provides information for
%     min   x'*[12 10;10 16]*x+[2 18]*x 
%     s.t. \|x\|_0<=1
% where  s=1 

switch fgH    
    case 'obj'  % objective function
    out = x'*[12 10;10 16]*x+[2 18]*x;     
    
    case 'grad'  % gradient 
    out = 2*[12 10;10 16]*x+[2; 18];    
    
    case 'hess'  % Hessian matrix 
    out = 2*[12 10;10 16];  
end
 
end


 
