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
from math import exp,log,sqrt,cos,pi,sin,pi,e
import matplotlib.pyplot as plt

# solution function
def func(x):
    y = (2 * exp(1))* x *(exp(-x)) - exp(x)
    # y = (25/(pi**2))*sin(pi*x) + x #example from class
    return y

# righthand side function
def rfunc(x):
    y = -(4 * exp(x))
    # y = 25 * sin(pi*x) #example from class
    return y

# coefficients of the differential equations
# A=-1; B=0; C=0 #example from class
A=1; B=2; C=1

# the domain [xa,xb]
# xa=0; xb=1 #example from class
xa=0; xb=2

# boudary conditions
alpha = func(xa); beta = func(xb)

# index for methods: 0:Thomas, 1:Jacobi, 2:Gauss-Sidel, 3: SOR
imethod=3

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
    u[0] = b[0]
    for k in range(1,n):
        l[k] = a[k]/u[k-1]
        u[k] = b[k] - l[k] * c[k-1]

    # Solve Lz = r
    z[0] = r[0]
    for k in range(1,n):
        z[k] = r[k] - l[k] * z[k-1]

    # Solve Uw = z.
    w[n-1] = z[n-1] / u[n-1]
    for k in range(n-2, -1, -1):
        w[k] = (z[k] - (c[k] * w[k+1]))/u[k]

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
    coA = (A/float(h**2)) - (B/float(2*h))
    coB = (float(-2*A)/float(h**2)) + C
    coC = (float(A)/float(h**2)) + (float(B)/float(2*h))

    # claim the vectors needed
    xh = zeros(n,float) #x-values
    yh = zeros(n,float) #true y-values at grids
    wh = zeros(n,float) #computed y-values at grids
    r=zeros(n,float)    #right handside of the equations

    #begin
    #assisgn values for xh,yh,r

    for i in range(0,n):
        xh[i] = xa+((i+1)*h)
        yh[i] = func(xh[i])

    #end


    r[0] = rfunc(xh[0]) - (coA*alpha)
    r[n-1]= rfunc(xh[n-1]) - (coC*beta)

    for i in range(1,n-1):
        r[i] = rfunc(xh[i])  # assign values in the middle


    # vectors needed for direct methods
    if (imethod == 0):
        # Thomas's algorithm
        a=zeros(n,float);b=zeros(n,float);c=zeros(n,float)
        # assign values for a,b,c
        for i in range(n):
            b[i] = coB
        for i in range(1, n):
            a[i] = coA
        for i in range(0, n-1):
            c[i] = coC

        wh = LU3315(a,b,c,r)

    else:
        # Iterative Methods
        tol=10**(-8); err=1 # initial error and error tolerance
        wh1=zeros(n,float) # vectors for computed values


        while (err>tol):
            icount[icase]=icount[icase]+1
            if (imethod==1):
                # Jacobi
                wh[0] = (r[0] - (coC * wh1[1]))/coB
                for i in range(1, n-1):
                    wh[i] = (r[i] - coC*wh1[i+1] - coA*wh1[i-1])/coB
                wh[n-1] = (r[n-1] - (coA * wh1[n-2])) / coB;

            elif (imethod==2):
                # Gauss-Seidel
                wh[0] = (r[0] - coC * wh[1]) / coB
                for i in range(1, n-1):
                    wh[i] = (r[i] - (coA * wh[i-1]) - (coC*wh[i+1]))/coB
                wh[n-1] = (r[n-1] - coA * wh[n-2]) / coB
            else:
                # SOR
                wh[0] = (coB*wh[0] + (wopt*r[0]) - (coC*wh[1]) - (coB*wh[0]))/coB
                for i in range(1,n-1):
                    wh[i] = (coB*wh[i] + wopt*(r[i] - coA*wh[i-1] - coC*wh[i+1] - coB*wh[i]))/coB
                wh[n-1] = (coB * wh[n-1] + (wopt * (r[n-1] - (coA * wh[n-2]) - (coB * wh[n-1])))) / coB

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
    kx=int(icase/2); ky=icase%2

    ax[kx,ky].plot(xe,ye,'-b',xplot,wplot,'ro' )
    ax[kx,ky].set_title('h='+str(h))


plt.savefig('result.pdf',format='pdf')
plt.show()
