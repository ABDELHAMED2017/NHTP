% demon compressed sensing problems with real data
clc; clear; close all;

load 'DrivFace.mat';
load 'nlab.mat';   %'identity.mat';

[m,n]     = size(A);
s         = 10;
data.A    = A;
data.b    = y;
data.At   = data.A';

pars.Draw = 1;
out       = NHTP(n,s,data,'CS',pars)  

fprintf('\nSample size:       m=%4d,n=%4d\n', m,n);
fprintf('CPU time:         %6.3fsec\n',  out.time);
fprintf('Objective value:   %5.3e\n\n', out.obj);

if isfield(pars,'Draw') && pars.Draw
saveas(figure(1), 'outputs\DriverFaces.fig');
saveas(figure(1), 'outputs\DriverFaces.eps','epsc');
end
