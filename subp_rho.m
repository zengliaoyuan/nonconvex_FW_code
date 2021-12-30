function [x, lambda] = subp_rho( w, xi,  sigmatilde, rho)
% This solves subproblems of FW for solving min ||Ax-b||^2/2  s.t. ||x||_1 + rho*||x||^2/2  - (mu*||x|| + rho*||x||^2/2)
% i.e, solves
% min w'*x
% s.t. ||x||_1+rho*||x||^2/2 - mu*||x^k||- rho*||x^k||^2/2 -  (mu*x^k/||x^k||+rho*x^k)'*(x-x^k) <= sigma
% input: w = A'*(Ax-b);  xi = mu*x^k/||x^k|| + rho*x^k;  sigmatilde = sigma - rho*||x^k||^2/2

n = length(w);

if norm(w) < 1e-8
    x = zeros(n, 1);
    return;
end

sigmahat = rho*sigmatilde;

[Case, sg] = CaseSg_rho(w, xi, n);    %computing each cases and segments of each case

p = sg(:);
[pp,I] = sort(p,'ascend');
[N,~] = size(I);
count = ones(n,1);
Fv = -0.1*ones(N,1);
theta = -0.1;     
a = sum(Case(:,1));
c = sum(Case(:,2));
ae =  sum(Case(:,5));
ce =  sum(Case(:,6));
Fv(N) = ae*(pp(N))^2 + ce - sigmahat;   % function value at the end break point
if Fv(N)<= 0
    theta = sqrt((sigmahat-ce)/ae);
else
    % computing function values on the break points
    for j=1:N
        Fv(j) = a*(pp(j))^2 + c - sigmahat;
        if Fv(j)>=0
            theta = sqrt((sigmahat-c)/a);
            break
        end
        k = mod(I(j),n);
        if k == 0
            k=n;
        end
        
        if count(k) <= 2
            a = a - Case(k,2*count(k) - 1) + Case(k,2*count(k) + 1);
            c = c - Case(k,2*count(k) ) + Case(k,2*count(k) + 2);
            count(k) = count(k) +1;
        end
    end
end

z = xi - w*theta;
xhat = max(0,abs(z)-1).*sign(z);
x = xhat/rho;
lambda = 1/theta;




