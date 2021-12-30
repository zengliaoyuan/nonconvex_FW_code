# nonconvex_FW_code

This package implements new nonconvex Frank-Wolfe (FW) type methods for minimizing optimization problems with a single **nonconvex constraint** as follows:

![problem](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D%20%5Cbegin%7Barray%7D%7Bcl%7D%5Cmin%20&%20f(x)%5C%5C%7B%5Crm%20s.t.%7D%20&%20P_1(x)%20-%20P_2(x)%5Cleq%20%5Csigma,%5Cend%7Barray%7D) 

<br />

## Matlab source codes are:


FW_cs: implementation of our nonconvex FW method and away-step method for problems with ![P_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D%20P_1) being convex <br />
fwstep: subroutine for the FW linear oracle in the FW method <br />
awstep: subrountine for the away-step linear oracle in the away-step FW method

<br />

FW_rho: implementation of our nonconvex FW method for problems with ![P_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D%20P_1) being strongly convex whose modulus is ![rho](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D%20%5Crho)  <br />
subp_rho: subroutine for the FW_rho code  <br />
CaseSg_rho: subroutine for the subp_rho code  <br />

<br />

runcodeFW: runcode. Can be run directly. Figures generated will be saved in "Results"  <br />
plotFvaltrace: subroutine for plotting in the runcodeFW code <br />
savefig: subroutine for saving figures in the plotFvaltrace code

<br />

## Reference

Please refer to ["Frank-Wolfe type methods for nonconvex
inequality-constrained problems"](https://arxiv.org/pdf/2112.14404.pdf) for more details.





