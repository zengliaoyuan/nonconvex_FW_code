function [Case, sg] = CaseSg_rho(w, xi, n)
%A subroutine of subp_rho

Case = zeros(n, 6);
sg = inf*ones(n, 2);
a_optn = w.^2/2;
BP1 = (xi + 1)./w;
BP2 = (xi - 1)./w;
c_optn1 = -(xi + 1).^2/2;
c_optn2 = -(xi - 1).^2/2;


% Case 1
I1 = (w == 0);
v1 = abs(xi(logical(I1))) - 1;
v2 = max(v1 , 0);
c1 = v2.*(v2/2-v1);
Case(logical(I1), [2, 4, 6]) = [c1,  c1, c1] ;

% Case 2
I2 = (w > 0).*(xi < -1);
a2 = a_optn(logical(I2));
c2 = c_optn1(logical(I2));
Case(logical(I2), :) = [a2 , c2, a2, c2, a2, c2];

% Case 3
I3 = (w > 0).*(xi > 1);
a3 = a_optn(logical(I3));
c31 = c_optn2(logical(I3));
c32 = c_optn1(logical(I3));
Case(logical(I3), [1, 2, 5, 6]) = [a3, c31, a3, c32];
sg(logical(I3), :) = [BP2(logical(I3)), BP1(logical(I3))];

% Case 4
I4 = (w > 0).*(abs(xi) <= 1);
a4 = a_optn(logical(I4));
c4 = c_optn1(logical(I4));
Case(logical(I4), [3, 4, 5, 6]) = [a4, c4, a4, c4];
sg(logical(I4), 1) = BP1(logical(I4));

% Case 5
I5 = (w < 0).*(xi < -1);
a5 = a_optn(logical(I5));
c51 = c_optn1(logical(I5));
c52 = c_optn2(logical(I5));
Case(logical(I5), [1, 2, 5, 6]) = [a5, c51, a5, c52];
sg(logical(I5), :) = [BP1(logical(I5)),  BP2(logical(I5))];

% Case 6
I6 = (w < 0).*(xi > 1);
a6 = a_optn(logical(I6));
c6 = c_optn2(logical(I6));
Case(logical(I6), :) = [a6, c6, a6, c6, a6, c6];

% Case 7
I7 = (w < 0).*(abs(xi) <= 1);
a7 = a_optn(logical(I7));
c7 = c_optn2(logical(I7));
Case(logical(I7), [3, 4, 5, 6]) = [a7, c7, a7, c7];
sg(logical(I7), 1) = BP2(logical(I7));

