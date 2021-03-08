from scipy.stats import norm
import scipy.integrate as integrate
import math

"Original code snipper"
# function payment = payoff(mu,sigma,p,pplus,pmin)
# %Inputs: mu and sigma of the agent (or coalition),
# %The p, pplus and pminus POU Tariff parameters
# %Output: Expected payment value
#
# if(pplus+pmin~=0)
#    ropt=pplus/(pplus+pmin);
# else
#     ropt=0;
# end
#
# fun = @(x) norminv(x);
# bigphi=0;
#
# if(ropt~=0)
# bigphi=integral(fun,0,ropt);
# end
#
# payment = p*mu - sigma*(pplus+pmin)*bigphi;
# end

def payoff(mu,sigma,p,pplus,pmin,x):
    "Added an x value for the @x, this will be declared outside"
    if pplus+pmin!=0:
        ropt = pplus / (pplus + pmin)
    else:
        ropt = 0
    fun = norm.ppf(x)
    bigphi = 0
    if ropt!=0:
        bigphi = integrate.quad(fun, 0, ropt)
    return p * mu - sigma * (pplus + pmin) * bigphi


"Original code snipper"
# function sd=sigma_set(vect)
# %Input: vector of individual sigmas
# %Output: sigma of coalition
# sz=size(vect);
# s=0;
# for(i=1:sz(2))
#   s=s+vect(i)^2;
# end
# sd=sqrt(s);
# end

def sd(vect):

    s= 0
    for i in range(0,len(vect[0])):
        s = s + vect(i) ^ 2

    return math.sqrt(s)

