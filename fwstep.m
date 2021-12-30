function ufw = fwstep(grad, xi, sigma)
%% Linear oracle. Solve:
% min <grad, u>
% s.t. \|u\|_1 - <xi, u> <= sigma 
% where xi = mu*Sgn(x)


%% main
Sgn =  min(2*sign(grad)+1, 1);      % sign of the gradient 
b = -abs(grad)./(1+xi.*Sgn);
[~, I] = min(b);                % min outputs 1 index even there is more than 1 maximal entries
n = size(xi, 1);
ufw = zeros(n, 1);
ai_sgn = Sgn(I);
ufw(I) = -sigma*ai_sgn /( 1 + xi(I)*ai_sgn );


