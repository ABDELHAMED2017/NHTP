function Out = NHTP(n, s, data,  probname, pars)

% This code aims at solving the sparsity constrained optimization with form
%
%         min_{x\in R^n} f(x)  s.t.  \|x\|_0<=s
%
% where s is the given sparsity, which is << n.  
%
% Inputs:
%     n       : Dimension of the solution x, (required)
%     s       : Sparsity level of x, an integer between 1 and n-1, (required)
%     data    : Input functions include the objective f(x), its gradient and Hessian
%               For CS and logistic regression, data also could be in the form:
%               data.A  --  m x n order measurement matrices 
%               data.At --  the transpose of data.A, i.e., data.At=data.A' 
%               data.y  --  m x 1 order observation vector 
%     probname: Name of problem, should be one of {'CS','SLR','SCO'}
%     pars:     Parameters are all OPTIONAL
%               pars.x0     --  Starting point of x,   pars.x0=zeros(n,1) (default)
%               pars.eta    --  A positive parameter,  a default one is given related to inputs  
%               pars.IterOn --  Results will  be shown for each iteration if pars.IterOn=1 (default)
%                               Results won't be shown for each iteration if pars.IterOn=0 
%               pars.draw   --  A  graph will be drawn if pars.draw=1 (default) 
%                               No graph will be drawn if pars.draw=0 
%
% Outputs:
%     Out.sol:           The sparse solution x
%     Out.sparsity:      Sparsity level of Out.sol
%     Out.normgrad:      L2 norm of the gradient at Out.sol  
%     Out.error:         Error used to terminate this solver 
%     Out.time           CPU time
%     Out.iter:          Number of iterations
%     Out.grad:          Gradient at Out.sol
%     Out.obj:           Objective function value at Out.sol 
%     Out.eta:           Final eta
%
%
%%%%%%%    Send your comments and suggestions to                     %%%%%%
%
%%%%%%%    shenglong.zhou@soton.ac.uk                                %%%%%%
% 
%%%%%%%    Warning: Accuracy may not be guaranteed!!!!!              %%%%%%
warning off;
t0     = tic;

if  nargin<3
    disp(' No enough inputs. No problems will be solverd!'); return;
elseif nargin==3
    disp(' We will solve this problem in a general way');
    funcNam = str2func('general_example');
elseif nargin>=4
    switch      probname
    case 'CS' ; funcNam = str2func('compressed_sensing');
    case 'SLR'; funcNam = str2func('logistic_regression');
    otherwise ; funcNam = str2func('general_example'); 
    end
end

% define the functions
funcs    = @(var,flag,T1,T2)funcNam(var,flag,T1,T2,data);
 
if nargin>=3
    if nargin<5; pars=[]; end
    if isfield(pars,'IterOn');iterOn = pars.IterOn; else; iterOn = 1;     end
    if isfield(pars,'Draw');  Draw   = pars.Draw;   else; Draw   = 1;     end
    if isfield(pars,'MaxIt'); ItMax  = pars.MaxIt;  else; ItMax  = 2000;  end
    if isfield(pars,'x0');    x0     = pars.x0;     else; x0 = zeros(n,1);end 
    
    [obj,g]  = funcs(x0,'ObjGrad',[],[]); 
    if isfield(pars,'eta')      
        eta  = pars.eta;       
    else
        [~,g1] = funcs(ones(n,1),'ObjGrad',[],[]) ;
        abg1   = abs(g1);
        T      = find(abg1>1e-8);
        maxe   = sum(1./(abg1(T)+eps))/nnz(T);
        if  isempty(T) 
            eta    = 10*(1+s/n)/min(10, log(n));
        else
            if maxe>2
            eta  = (log2(1+maxe)/log2(maxe))*exp((s/n)^(1/3));
            elseif maxe<1
            eta  = (log2(1+ maxe))*(n/s)^(1/2);    
            else
            eta  = (log2(1+ maxe))*exp((s/n)^(1/3));
            end     
        end
    end      
end
                 
x       = x0;
beta    = 0.5;
tol     = 1e-6; 
sigma   = 5e-5;
I       = 1:n;
delta   = 1e-10;
pcgtol  = tol*s;
pcgs0   = 2000;
T0      = zeros(s,1);
Error   = zeros(1,ItMax);
Obj     = zeros(1,ItMax);

fprintf(' Start to run the sover...\n'); 
if iterOn 
fprintf('\n Iter             Error                Ojective \n'); 
fprintf('--------------------------------------------------------\n');
end

% Initial check for the starting point
if FNorm(g)==0 
fprintf('Starting point is a good stationary point, try another one...\n'); 
return;  
end


if max(isnan(g))
x0      = zeros(n,1);
rind    = randi(n);
x0(rind)= rand;
[obj,g] = funcs(x0,'ObjGrad',[],[]);
end
 
% The main body  
for iter = 1:ItMax
     
    xtg   = x0-eta*g;
    [~,T] = maxk(abs(xtg),s); 
    flag  = isempty(setdiff(T,T0));          
    Tc    = setdiff(I,T);
    
    % Calculate the error for stopping criteria 
    
    xtaus       = max(0,max(abs(g(Tc)))-min(abs(x(T)))/eta);
    if flag
    FxT         = sqrt(FNorm(g(T)));
    Error(iter) = xtaus + FxT;
    else
    FxT         = sqrt(FNorm(g(T))+ FNorm(x(Tc)));
    Error(iter) = xtaus + FxT;    
    end
     
    if iterOn 
    fprintf('%4d             %5.2e              %5.2e\n',iter,Error(iter),obj); 
    end
             
    % Stopping criteria
    if Error(iter)<tol; break;  end
    
    
    % update next iterate
    if  iter   == 1 || flag           % update next iterate if T==supp(x^k)     
        H      =  funcs(x0,'Hess',T,[]); 
        if s   < pcgs0
        d      = -H\g(T);
        else
        d      = -pcg(H,g(T),pcgtol,50); 
        end
        dg     =  sum(d.*g(T));
        if dg  >  max(-delta*FNorm(d), -FNorm(g(T)))
        d      = -g(T); 
        end
    else                              % update next iterate if T~=supp(x^k)
        TTc    = intersect(T0,Tc); 
        [H,D]  = funcs(x0,'Hess',T,TTc);
        if s   < pcgs0
        d      = H\( D*x0(TTc)- g(T));
        else
        d      = pcg(H,D*x0(TTc)- g(T), pcgtol,50); 
        end
        Fnz    = FNorm(x(TTc))/4/eta;
        dgT    = sum(d.*g(T));
        dg     = dgT-sum(x0(TTc).*g(TTc));
        
        delta0  = delta;
        if Fnz > 1e-4; delta0 = 1e-4; end
 
        if dgT > max(-delta0*FNorm(d)+Fnz, -FNorm(g(T)))
        d      = -g(T);  
        end        
    end
    

    alpha    = 1; 
    x(Tc)    = 0;    
    obj0     = obj;        
    Obj(iter)= obj;
    
    % Amijio line search
    for i=1:5
        x(T)   = x0(T) + alpha*d;
        obj    = funcs(x,'ObjGrad',[],[]);
        if obj < obj0+alpha*sigma*dg; break; end        
        alpha  = beta*alpha;
    end
    
    % Hard Thresholding Pursuit if the obj increases
    if obj  > obj0 
       x(T) = xtg(T); 
       obj  = funcs(x,'ObjGrad',[],[]); 
    end

    % Stopping criteria
    if  iter>10 && abs(obj - obj0)<1e-6*(1+abs(obj)) 
        if obj > obj0
           iter    = iter-1; 
           x       = x0; 
           T       = T0; 
           [obj,g] = funcs(x,'ObjGrad',[],[]);
        end   
        break;
     end 
 
    T0      = T; 
    x0      = x; 
    [obj,g] = funcs(x,'ObjGrad',[],[]);
    
    % Update tau
    if  mod(iter,50)==0  
        if Error(iter)>1/iter^2  
        if iter<1500; eta = eta/1.05; 
        else;         eta = eta/1.5; 
        end     
        else;         eta = eta*1.05;   
        end
    end     
end
x(abs(x)<1e-4)=0;

% results output
time        = toc(t0);
Out.sparsity= nnz(x);
Out.normgrad= sqrt(FNorm(g)); 
Out.error   = sqrt(FNorm(g(T))+ FNorm(x(setdiff(I,T))));
Out.time    = time;
Out.iter    = iter;
Out.grad    = g;
Out.sol     = x;
Out.obj     = obj; 
Out.eta     = eta;

if Draw
    figure
    subplot(121) 
    Obj(iter)= obj;  
    PlotFun(Obj,iter,'r.-','f(x^k)'); 
    subplot(122) 
    PlotFun(Error,iter,'r.-','error') 
end

fprintf('\n--------------------------------------------------------\n');

if  Out.normgrad<1e-5
fprintf(' A global optimal solution might be found\n');
fprintf(' because of ||gradient||=%5.2e!\n',Out.normgrad); 
if Out.iter>1500
fprintf('\n Since the number of iterations reaches to %d\n',Out.iter);
fprintf(' Try to rerun the solver with a smaller pars.eta=%5.2e\n',Out.eta); 
end
fprintf('--------------------------------------------------------\n');
end

end

function Fnorm = FNorm(x)
Fnorm = sum(x.*x);
end


function  maxeta = GetInitEta(funcs,x0)
[~,g1]  = funcs(x0,'ObjGrad',[],[]);
 abg1   = abs(g1);
T       = find(abg1>1e-8);
if isempty(T)
maxeta  = 0;
else
maxeta  = sum(abs(x0(T))./(abg1(T)+eps))/nnz(T);
end
end