# nonconvex_FW_code

This package implements new nonconvex Frank-Wolfe (FW) type methods for minimizing optimization problems with a single **nonconvex constraint** as follows:

![problem](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Bequation*%7D%0A%5Cbegin%7Barray%7D%7Bcl%7D%0A%5Cmin%20%26%20f(x)%20%5C%5C%0A%7B%5Crm%20s.t.%7D%20%26%20P_1(x)%20-%20P_2(x)%20%5Cleq%20%5Csigma%0A%5Cend%7Barray%7D%0A%5Cend%7Bequation*%7D)



## Matlab source codes are:


FW_cs: implementation of our nonconvex FW method and away-step method for problems with ![P_1](https://render.githubusercontent.com/render/math?math=P_1) being convex <br />
fwstep: subroutine for the FW linear oracle in the FW method <br />
awstep: subrountine for the away-step linear oracle in the away-step FW method

<br />

FW_rho: implementation of our nonconvex FW method for problems with ![P_1](https://render.githubusercontent.com/render/math?math=P_1) being strongly convex whose modulus is ![rho](https://render.githubusercontent.com/render/math?math=\rho)  <br />
subp_rho: subroutine for the FW_rho code  <br />
CaseSg_rho: subroutine for the subp_rho code  <br />

<br />

runcodeFW: runcode. Can be run directly. Figures generated will be saved in "Results"  <br />
plotFvaltrace: subroutine for plotting in the runcodeFW code <br />
savefig: subroutine for saving figures in the plotFvaltrace code


## Reference

Please refer to ["Frank-Wolfe type methods for nonconvex
inequality-constrained problems"](https://arxiv.org/pdf/2112.14404.pdf) for more details.





