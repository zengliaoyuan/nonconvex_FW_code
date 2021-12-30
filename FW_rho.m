function [x, iter, fval, trace] = FW_rho(A, b, mu, sigma, opts)

%% This uses FW to solve problems with P_1 strongly convex:
% min  || Ax - b ||^2/2
% s.t.  ||x||_1+rho*||x||^2/2 - (mu||x||+rho*||x||^2/2) <= sigma

%% % inputs
%  A  --  the sensing matrix
%  b  --  the observations
%  sigma  --  positive number
%  opts  -- options for the parameters used in the algorithms with fields below

%  fields of opts
%       opts.x0  --  a feasible initial point for the algorithm; The default value is zeros(n, 1).
%       opts.tol  -- the tolerance for termination; Default value is 1e-2.
%       opts.maxiter -- maximal iteration number; Default value is inf.
%       opts.isbdryb  -- nonnegative integer; frequence for boundary boosting;  if 0, then no boundary boosting; default isbdryb = 1;
%       opts.disp  -- nonnegative integer; frequence for displaying information; if 0, no display; default disp = 0
%       opts.maxtime -- the maximun allowed cputime; default value is inf
%       opts.rho  -- the module used for constructing the strongly convex P_1; default rho = 1

%% % Initialize
tstart = tic;
n = size(A, 2);

if isfield(opts, 'rho')
    rho = opts.rho;
else
    rho = 1;
    warning('The strongly convex modulus is 1 by default. \n')
end

if isfield(opts, 'isbdryb')
    isbdryb = opts.isbdryb;
else
    isbdryb = 1;
end

if isfield(opts, 'maxtime')
    maxtime = opts.maxtime;
else
    maxtime = inf;
end

if isfield(opts, 'maxiter')
    maxiter = opts.maxiter;
else
    maxiter = inf;
end

if isfield(opts, 'tol')
    tol = opts.tol;
else
    tol = 1e-2;
end

if isfield(opts, 'disp')
    disp = opts.disp;
else
    disp = 0;
end

if isfield(opts, 'x0')
    x = opts.x0;
else
    x = zeros(n, 1);        % zeros(n, 1) is feasible to the constraint ||x||_1 - mu||x|| <= sigma
end

c = 1e-4;
alpha_init = 1;
iter = 1;
traceTime = [];
traceFval = [];

Axb = A*x - b;
grad = A'*Axb;
fval = norm(Axb)^2/2;

fprintf(' ****************** Start of FW_rho with rho = %g********************\n', rho)
fprintf('fval = %6.6e at the initial point \n', fval);
fprintf('  iter     iter1      fval           fval1           c_res             fdiff             alpha      dderiv \n')

%% % main loop
while 1==1
    if norm(x) < 1e-8
        xi = zeros(n, 1);
    else
        xi = (mu/norm(x) + rho)*x;
    end
    sigmatilde = sigma - norm(x)^2*(rho/2);
    
    %solve the subproblem    
    u = subp_rho(grad, xi, sigmatilde, rho);
    dir = u - x;
    dderiv = grad'*dir;
    
    timeNow = toc(tstart);
    traceTime = [traceTime; timeNow];
    traceFval = [traceFval; fval];
    
    %% % termination
    if abs(dderiv)/max(abs(fval), 1) < tol || iter > maxiter  || toc(tstart) > maxtime
        fprintf(' Termination of FW_rho:  iter =%d, nnz = %d, fval =%g, c_res = %6.4e, dderiv = %6.4e \n', iter, nnz(abs(x) > 1e-6), fval, (norm(x,1) - mu*norm(x) - sigma)/sigma, dderiv)
        break
    end
    
    %% % linesearch
    alpha = alpha_init;
    iter_ls = 0;
    Adir = A*dir;
    while 1 == 1
        xtilde = x + alpha*dir;
        Axbtilde = Axb + alpha*Adir;
        ftilde = norm(Axbtilde)^2/2;
        if ftilde <= fval + c*alpha*dderiv
            break;
        else
            alpha = alpha/2;
            iter_ls = iter_ls +1;
        end
    end
    
    c1 = norm(xtilde,1) - mu*norm(xtilde);
    c_res = (c1-sigma)/sigma;
    Fdiff = abs(ftilde - fval)/ max(1, fval);
    if disp && mod(iter,disp) == 0
        fprintf(' %5d|%5d| %6.6e|%6.6e|%6.3e|%6.3e|%6.3e|%6.3e\n',iter, iter_ls, fval, ftilde, c_res, Fdiff, alpha, dderiv)
    end
   
    
    %% % boundary boosting  and updation
    if isbdryb && mod(iter, 1) == 0 && c1 > 1e-8
        ratio = sigma/c1;
        x = ratio*xtilde;
        Axb =  ratio*(Axbtilde + b) - b;
        fval = norm(Axb)^2/2;
        if fval > ftilde
            x = xtilde;
            Axb = Axbtilde;
            grad = A'*Axb;
            fval = ftilde;
        else
            grad = A'*Axb;
        end
    else
        x = xtilde;
        Axb = Axbtilde;
        grad = A'*Axb;
        fval = ftilde;
    end
    
    iter = iter +1;
        
    if iter_ls >= 1
        alpha_init = min(max(alpha, 1e-8), 1);
    else
        alpha_init = min(max(alpha*2, 1e-8), 1);
    end
end

trace = [traceTime, traceFval];

end


