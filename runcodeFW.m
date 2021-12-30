clear;
clc;

rand('seed', 2012);
randn('seed', 2011);

opts.tol = 1e-10;
opts.maxiter = 10000;
opts.disp = 1000;   % the frequence of displaying iterative information; 0: no display; >0: display every disp iterates
opts.isbdryb = 1;      % use the boundary boosting technique for all methods


repeat = 30;
Aindex = [1 3];
maxtime_arr = [5 30];
Arho = [0.01  0.05   0.1];

trace_all = cell(length(Aindex), 1);
for ii = 1:length(Aindex)
    index = Aindex(ii);
    opts.maxtime = maxtime_arr(ii);
    
    m = 720*index;
    n = 2560*index;
    k = 80*index;  % originally 80
    opts.x0 = zeros(n,1);       %  use the same initial point for all methods
    
    trace_index = cell(repeat, 1);
    for rr = 1:repeat
        fprintf('\n\n index = %d, repeat = %d \n', index, rr)
        trace_index{rr} = cell(2+length(Arho), 1);
        A = randn(m, n);
        for j = 1 : n
            A(:, j) = A(:, j)/norm(A(:, j));
        end
        
        % generate the original signal
        I = randperm(n);
        J = I(1 : k);
        xorig = zeros(n, 1);
        xorig(J) = randn(k, 1);
        
        noise = 0.01*randn(m, 1);
        b = A*xorig + noise;
        
        mu = 0.5;
        sigma = (norm(xorig,1)  - mu*norm(xorig))/1.1;
        
        %% % test the Frank-Wolfe method with boundary boosting every iteration
        opts.isafw = 0;
        tstart0 = tic;
        [x_fw, iter_fw, fval_fw, trace_fw] = FW_cs(A, b, mu, sigma, opts);
        cputime_fw = toc(tstart0);
        
        trace_index{rr}{1} = trace_fw;
        
        RecErr_fw = norm(x_fw - xorig)/norm(xorig);
        Residual_fw = (norm(x_fw, 1) - mu*norm(x_fw) - sigma)/sigma; % relative feasibility
        nnz_fw = nnz(abs(x_fw) > 1e-6);
        fprintf(' FW with boundary operation:  time = %6.4f, iter =%d, nnz = %g, fval = %g, c_res = %g, RecErr = %g\n',...
            cputime_fw, iter_fw, nnz_fw, fval_fw, Residual_fw, RecErr_fw)
       
        
        %% % test the Away-Step FW method with boundary boosting every iteration
        opts.isafw = 1;
        tstart0 = tic;
        [x_afw, iter_afw, fval_afw, trace_afw] = FW_cs(A, b, mu, sigma, opts);
        cputime_afw = toc(tstart0);
        
        trace_index{rr}{2} = trace_afw;
        
        RecErr_afw = norm(x_afw - xorig)/norm(xorig);
        Residual_afw = (norm(x_afw, 1) - mu*norm(x_afw) - sigma)/sigma;
        nnz_afw = nnz(abs(x_afw) > 1e-6);
        fprintf(' AFW with boundary operation: time = %6.4f, iter =%d, nnz = %g, fval = %g, c_res = %g, RecErr = %g\n',...
            cputime_afw, iter_afw, nnz_afw, fval_afw, Residual_afw, RecErr_afw)
        
        
        %%  % test the FW_rho method
        for kk = 1:length(Arho)
            opts.rho = Arho(kk);
            tstart = tic;
            [x_rho, iter_rho, fval_rho, trace_rho] = FW_rho(A, b, mu, sigma,opts);
            cputime_rho = toc(tstart);
            
            trace_index{rr}{2+kk} = trace_rho;
            
            RecErr_rho = norm(x_rho - xorig)/norm(xorig);
            Residual_rho = (norm(x_rho, 1) - mu*norm(x_rho) - sigma)/sigma;
            nnz_rho= nnz(abs(x_rho) > 1e-6);
            fprintf(' FW_rho with rho %g: time = %6.4g, iter =%d, nnz = %g, fval = %g, c_res = %g, RecErr = %g \n',...
                opts.rho, cputime_rho, iter_rho, nnz_rho, fval_rho, Residual_rho, RecErr_rho)
        end
    end
    trace_all{ii} = trace_index;
end

a = clock;
matfilename = ['data_cs' '-' date '-' int2str(a(4)) '-' int2str(a(5)) '.mat'];
save(matfilename,  'trace_all',   'Aindex', 'maxtime_arr',  'repeat',  'Arho');
save haha trace_all Aindex maxtime_arr repeat Arho

plotFvaltrace;





