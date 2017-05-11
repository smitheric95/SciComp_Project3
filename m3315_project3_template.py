"""
This function solves a two-point boundary value problem
Ay''+By'+Cy=r, y(a)=alpha, y(b)=beta,
using a finite-difference scheme.
The discreitzed system is solved by direct tridiagonal LU method
and iterative methods including Jacobi, Gauss-Seidel, and SOR
Note:
1. Replace single line *** by single line code and
relace double line *** by multiple lines of code.
2. Python array index start with 0.
"""

import numpy
from numpy import exp,zeros,arange,absolute,max
from math import exp,log,sqrt,cos,pi,sin,pi
import matplotlib.pyplot as plt

# solution function
def func(x):
    y=***
    return y

# righthand side function
def rfunc(x):
    y=***
    return y

# coefficients of the differential equations
A=***; B=***; C=***

# the domain [xa,xb]
xa=***; xb=***

# boudary conditions
alpha = func(***); beta = func(***)

# index for methods: 0:Thomas, 1:Jacobi, 2:Gauss-Sidel, 3: SOR
imethod=0

# number of cases, each case has different # of unknowns
ncase=4

# table for results
tbErr=zeros((ncase,4),float)

# the counter for iterative methods
icount=zeros(ncase,int)

#LU factorization based on Thomas's algorithm
#input: a, b, c, r - matrix elements and right hand side vector
#output: w - solution of linear system
def LU3315(a,b,c,r):
    n = len(r)
    w = zeros(n,float)
    l = zeros(n,float)
    u = zeros(n,float)
    z = zeros(n,float)
    #Determine L,U factors
    ***
    ***

    # Solve Lz = r
    ***
    ***

    # Solve Uw = z.
    ***
    ***

    return w


# the main code starts here
fig, ax = plt.subplots(2,2,sharex=True, sharey=True)

for icase in range(ncase):
    n = 2**(icase+2)-1   #number of unknowns
    h = (xb-xa)/(n+1);  #mesh size
    wopt = 2/(1+sqrt(1-cos(pi*h)**2)) #optimized omega

    # exact value at a fine mesh
    d=0.0025
    xe = arange(xa,xb+d,d)
    ye = xe.copy()
    for i in range(len(xe)):
        ye[i] = func(xe[i])


    # matrix entry on tri-dialgonals
    coA=***
    coB=***
    coC=***


    # claim the vectors needed
    xh = zeros(n,float) #x-values
    yh = zeros(n,float) #true y-values at grids
    wh = zeros(n,float) #computed y-values at grids
    r=zeros(n,float)    #right handside of the equations

    #begin
    #assisgn values for xh,yh,r
    ***
    ***
    #end

    # vectors needed for direct methods
    if (imethod == 0):
        # Thomas's algorithm
        a=zeros(n,float);b=zeros(n,float);c=zeros(n,float)
        # assign values for a,b,c
        ***
        ***
        wh = LU3315(a,b,c,r)
    else:
        # Iterative Methods
        tol=10**(-8); err=1 # initial error and error tolerance
        wh1=zeros(n,float) # vectors for computed values


        while (err>tol):
            icount[icase]=icount[icase]+1
            if (imethod==1):
                # Jacobi
                ***
                ***

            elif (imethod==2):
                # Gauss-Seidel
                ***
                ***
            else:
                # SOR
                ***
                ***

            err=max(absolute(wh1-wh))
            wh1=wh.copy()

    #output
    tbErr[icase,0] = h
    tbErr[icase,1] = max(absolute(yh-wh))

    if (icase > 0):
        tbErr[icase,2] = (tbErr[icase-1,1]/tbErr[icase,1])
        tbErr[icase,3] = log(tbErr[icase,2],2)

    if (imethod != 0):
        print('case ',icase,'iteration number= ',icount[icase])

    if (icase==ncase-1):
        print(tbErr)

    # plot: you don't need to change anything here
    xplot = zeros(n+2,float); wplot = zeros(n+2,float)
    xplot[0]=xa; xplot[n+1]=xb
    wplot[0]=alpha; wplot[n+1]=beta
    for i in range(1,n+1):
        xplot[i] = xh[i-1]
        wplot[i] = wh[i-1]
    kx=icase/2; ky=icase%2
    ax[kx,ky].plot(xe,ye,'-b',xplot,wplot,'ro' )
    ax[kx,ky].set_title('h='+str(h))

plt.savefig('result.pdf',format='pdf')
plt.show()
