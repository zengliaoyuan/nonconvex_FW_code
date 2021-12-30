function  [x, iter, fval, trace] = FW_cs(A, b, mu, sigma, opts)
%% % This uses nonconvex Frank-Wolfe method to solve
% %   the optimization problem arising from compressed sensing:
% %   min_{x\in\R^n} ||Ax-b||^2/2  s.t. ||x||_1 - mu*||x|| <= sigma, where 0< mu <1

%% % inputs
%  A  --  the sensing matrix
%  b  --  the observations
%  sigma  --  positive number
%  opts  -- options for the parameters used in the algorithms with fields below

%  fields of opts
%       opts.x0  --  a feasible initial point for the algorithm; The default value is zeros(n, 1).
%       opts.tol  -- the tolerance for termination; Default value is 1e-2.
%       opts.maxiter -- maximal iteration number; Default value is inf.
%       opts.isbdryb  -- frequence for boundary boosting;  nonnegative integer; if 0, then no boundary boosting; default isbdryb = 1;
%       opts.isafw  -- 0/1; tag for the away-step operation; default value is 1
%       opts.disp  -- frequence for displaying information;  nonnegative integer; if 0, no display; default disp = 0
%       opts.maxtime -- the maximum allowed cputime; default value is inf

%% %  Initialization
tstart = tic;

n = size(A, 2);

if isfield(opts, 'maxtime')
    maxtime = opts.maxtime;
else
    maxtime = inf;
end

if isfield(opts, 'x0')
    x = opts.x0;
else
    x = zeros(n, 1);        % zeros(n, 1) is feasible to the constraint ||x||_1 - mu||x|| <= sigma
end

if isfield(opts, 'tol')
    tol = opts.tol;
else
    tol = 1e-2;
end

if isfield(opts, 'maxiter')
    maxiter = opts.maxiter;
else
    maxiter = inf;
end

if isfield(opts, 'disp')
    disp = opts.disp;
else
    disp = 0;
end

if isfield(opts, 'isbdryb')
    isbdryb = opts.isbdryb;
else
    isbdryb = 1;
end

if isfield(opts, 'isafw')
    isafw = opts.isafw;
else
    isafw = 1;
end

c = 1e-4;
epsilon = 1e-5;
zeta = 1e5;
iter = 1;
alpha_init = 1;
alpha_fw = alpha_init;
traceTime = [];
traceFval = [];

Axb = A*x - b;
grad = A'*Axb;
fval = norm(Axb)^2/2;

fprintf('fval = %6.6e at the initial point \n', fval)
fprintf('  iter     iter1      fval           fval1           c_res             fdiff             alpha      dderiv       is_aw  \n')

%% % main loop
while 1
    % find a subgradient of the constraint function
    if norm(x) < 1e-8
        xi = zeros(n,1);
    else
        xi = (mu/norm(x))*x;
    end
    
    % find the frank-wolfe direction
    u_fw = fwstep(grad, xi, sigma);
    d_fw = u_fw - x;
    dderiv_fw = grad'*d_fw;
    
    timeNow = toc(tstart);
    traceTime = [traceTime; timeNow];
    traceFval = [traceFval; fval];
        
    % termination
    if abs(dderiv_fw)/max(abs(fval), 1) < tol || iter > maxiter || toc(tstart) > maxtime
        if ~isafw
            fprintf(' Termination of FW: iter =%d, nnz = %g, fval =%g, cst_res = %g, dderiv = %6.4e \n', iter, nnz(abs(x) > 1e-6), fval, (norm(x,1) - mu*norm(x) - sigma)/sigma, dderiv);
        else
            fprintf(' Termination of AFW: iter =%d, nnz = %g, fval =%g, cst_res = %g,  dderiv = %6.4e \n', iter, nnz(abs(x) > 1e-6), fval, (norm(x,1) - mu*norm(x) - sigma)/sigma, dderiv);
        end
        break;
    end
    
    
    if isafw
        % find the away-step direction
        if norm(x) < 1e-8
            d_aw = zeros(n, 1);
            dderiv_aw = 0;
            alpha_aw = 0;
        else
            [d_aw, alpha_aw] = awstep(grad, x, sigma, mu, zeta);
            dderiv_aw = grad'*d_aw;
        end
        
        % choosing a direction & set alpha_init
        if dderiv_fw > dderiv_aw  && alpha_aw > epsilon
            is_aw = 1;
            dir = d_aw;
            dderiv = dderiv_aw;
            alpha_init = alpha_aw;
        else
            is_aw = 0;
            dir  = d_fw;
            dderiv = dderiv_fw;
        end
    else
        is_aw = 0;
        dir = d_fw;
        dderiv = dderiv_fw;
    end
    
    % line search
    iter_ls = 0;
    alpha = alpha_init;
    Adir = A*dir;
    while 1 == 1
        xtilde = x + alpha*dir;
        Axbtilde = Axb + alpha*Adir;
        ftilde = norm(Axbtilde)^2/2;
        
        if ftilde <= fval + c*alpha*dderiv
            break;
        else
            alpha = alpha/2;
            iter_ls = iter_ls + 1;
        end
    end
    
    c1 = norm(xtilde,1) - mu*norm(xtilde);
    c_res = (c1-sigma)/sigma;
    Fdiff = abs(fval - ftilde)/max(1, abs(fval));  
    if disp && mod(iter,disp) == 0
        fprintf(' %5d|%5d| %6.6e|%6.6e|%6.3e|%6.3e|%6.3e|%6.3e |%3d\n', iter, iter_ls, fval, ftilde, c_res, Fdiff, alpha, dderiv, is_aw)
    end

    %% %  boundary boosting and updation
    if isbdryb && mod(iter, 1) == 0  && c1 > 1e-8
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
    
    iter = iter + 1;
    
    if ~is_aw
        alpha_fw = alpha;
    end
    if (~is_aw) && iter_ls == 0 
        alpha_init = min(max(alpha_fw*2, 1e-8), 1);
    else
        alpha_init = min(max(alpha_fw, 1e-8), 1);
    end
end

trace = [traceTime, traceFval];

end


