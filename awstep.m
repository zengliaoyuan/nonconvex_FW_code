function [d_aw, alpha_aw] = awstep(a, x, sigma, mu, zeta)
%% Away-step oracle. Solve
% max <a, u>
% s.t. u \in S(x, xi)
% where xi = mu*Sgn(x)


%% main
I = find( x~=0 );
a_I = a(I);
x_I = x(I);
xi_I = (mu/norm(x)) * x_I;
u_I = sigma * (sign(xi_I) ./ (1 - abs(xi_I)));     %  note that xi_I\neq 0; so it is fine to use the sign function here
c_I = x_I ./ u_I;

b = a_I .* u_I;
[bmax, id_I] = max(b);

u_id = u_I(id_I);
x_id = x_I(id_I);
id = I(id_I);
d_aw = x;
d_aw(id) = x_id - u_id;
c_id = c_I(id_I);
alpha_aw = c_id / (1 - c_id);    

if sum(c_I) < 1
    xi_id = xi_I(id_I);
    uprime = - sigma * sign(xi_id) / (1 + abs(xi_id));          % choose the opposite direction of u_id
    cprime = (1 - sum(c_I)) * (1 + abs(xi_id) )/ 2;                 % compute the coefficient of the additional node uprime
    c_id = c_id + (1 - sum(c_I)) * (1 - abs(xi_id) ) / 2;           % update the coefficient of u_id
    bprime = a_I(id_I) * uprime;
    if bmax <= bprime
        d_aw = x;
        d_aw(id) = x_id - uprime;
        alpha_aw = cprime / (1 - cprime);
    else
        alpha_aw = c_id / (1 - c_id);
    end
end

alpha_aw = min(alpha_aw, zeta);

end