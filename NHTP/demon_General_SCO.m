% demon a general sparsity constrained problem
%     min    x'*[6 5;5 8]*x+[1 9]*x-sqrt(x'*x+1)  
%     s.t. \|x\|_0<=s
% where s=1
% 
clc; clear; close all;

n    = 2;    
s    = 1;

% you can find this function in 'examples'-->'general_sco'
data = @(var,flag)simple_ex4_func(var,flag);  
pars.Draw   = 1;
out         = NHTP(n,s,data,'SCO',pars)  
fprintf('\nProblem dimension: n=%d\n', n);
fprintf('CPU time:         %6.3fsec\n',  out.time);
fprintf('Objective value:  %5.2e\n\n', out.obj);

if isfield(pars,'Draw') && pars.Draw
saveas(figure(1), 'outputs\GenSCO.eps','epsc');
saveas(figure(1), 'outputs\GenSCO.fig');
end

 