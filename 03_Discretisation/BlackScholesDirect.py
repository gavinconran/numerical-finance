from math import *

# Cumulative normal distribution function
def CND(X):
    (a1,a2,a3,a4,a5) = (0.31938153, -0.356563782, 1.781477937, -1.821255978, 1.330274429)
    L = abs(X)
    K = 1.0/(1.0+0.2316419*L)
    w = 1.0-1.0/sqrt(2.0*pi)*exp(-L*L/2.0)*(a1*K+a2*K*K+a3*pow(K,3)+a4*pow(K,4)+a5*pow(K,5))
    if X<0:
        w = 1.0-w
    return w

def phi(x):
    return exp(-x*x/2.)/sqrt(2*pi)
   
    
#### Black Sholes Function
def BlackScholesCore(CallPutFlag,DF,F,X,T,v):
    vsqrt=v*sqrt(T)
    d1 = (log(F/X)+(vsqrt*vsqrt/2.))/vsqrt
    d2 = d1-vsqrt
    if CallPutFlag:
        return DF*(F*CND(d1)-X*CND(d2))
    else:
        return DF*(X*CND(-d2)-F*CND(-d1))
        
#### Black Sholes Vega
def BlackScholesVegaCore(DF,F,X,T,v):
	vsqrt=v*sqrt(T)
	d1 = (log(F/X)+(vsqrt*vsqrt/2.))/vsqrt
	return F*phi(d1)*sqrt(T)/DF
	
def BlackScholes(CallPutFlag,S,X,T,r,d,v):
    return BlackScholesCore(CallPutFlag,exp(-r*T),exp((r-d)*T)*S,X,T,v)
        
        
def volimpl(CallPutFlag,S,X,T,r,d,price,tolerance=1e-10):
    upperPrice=BlackScholes(CallPutFlag,S,X,T,r,d,1.0)
    lowerPrice=BlackScholes(CallPutFlag,S,X,T,r,d,0.0001)
    upperVol,lowerVol=5.0,0.0
    if price>upperPrice or price<lowerPrice:
        return -1
    while upperVol>lowerVol+tolerance:
        midVol=0.5*(upperVol+lowerVol)
        midPrice=BlackScholes(CallPutFlag,S,X,T,r,d,midVol)
        if midPrice<price:
            lowerVol=midVol
        else:
            upperVol=midVol
    return upperVol

def volimplCore(CallPutFlag,DF,F,X,T,price,tolerance=1e-4):
    upperPrice=BlackScholesCore(CallPutFlag,DF,F,X,T,1.0)
    lowerPrice=BlackScholesCore(CallPutFlag,DF,F,X,T,0.0001)
    upperVol,lowerVol=1.0,0.0
    if price>upperPrice or price<lowerPrice:
        return -1
    while upperVol>lowerVol+tolerance:
        midVol=0.5*(upperVol+lowerVol)
        midPrice=BlackScholesCore(CallPutFlag,DF,F,X,T,midVol)
        if midPrice<price:
            lowerVol=midVol
        else:
            upperVol=midVol
    return upperVol

if __name__=="__main__":
	S=100.
	X=100.
	r=0.
	v=0.2
	d=0.
	T=1.
	DF=exp(-r*T)
	F=exp((r-d)*T)*S
	vsqrt=v*sqrt(T)
	d1 = (log(F/X)+(vsqrt*vsqrt/2.))/vsqrt
	d2 = d1-vsqrt
	print BlackScholes(True,S,X,T,r,d,v)

