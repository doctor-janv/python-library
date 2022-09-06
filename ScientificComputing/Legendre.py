import numpy as np
import math


#============================
def Legendre(N,x):
    Pnm1 = 1
    Pn   = x
    Pnp1 = 0
    
    if (N==0):
        return 1
    
    if (N==1):
        return x
    
    for n in range(2,N+1):
        ns=n-1
        Pnp1 = ((2*ns+1)/(ns+1))*x*Pn -(ns/(ns+1))*Pnm1
        Pnm1 = Pn
        Pn = Pnp1
        
    return Pnp1  

#============================
def dLegendredx(N,x):
    if (N==0):
        return 0
    
    if (N==1):
        return 1
    
    return (N*x/(math.pow(x,2)-1))*Legendre(N,x)- \
           (N/(math.pow(x,2)-1))*Legendre(N-1,x)
           
#============================
def LegendreRoots(N,maxiters=1000,tol=1.0e-10):
    xn = np.linspace(-0.999,0.999,N)  #Initial guessed values

    wn = np.zeros((N))
    
    
    for k in range(0,N):
        i=0
        while (i<maxiters):
            xold = xn[k]
            a = Legendre(N,xold)
            b = dLegendredx(N,xold)
            c = 0
            for j in range(0,k):
                c=c+(1/(xold-xn[j]))
            
            xnew = xold - (a/(b-a*c))
            
            res=abs(xnew-xold)
            xn[k]=xnew

            if (res<tol):
                break
            i=i+1
        
        wn[k] = 2*(1-xn[k]*xn[k])/(N+1)/(N+1)/ \
                Legendre(N+1,xn[k])/Legendre(N+1,xn[k])
    
    #======================Sorting the roots
    for i in range(0,N-1):
        for j in range(0,N-i-1):
            if (xn[j] > xn[j+1]):
                tempx = xn[j+1]
                tempw = wn[j+1]
                xn[j+1]=xn[j]
                wn[j+1]=wn[j]
                xn[j]=tempx
                wn[j]=tempw
                
    #for k in range(0,N):
    #    print("Finding root %d of %d" % (k,N), end='')
    #    print(" root %f, weight=%f, test=%f" %(xn[k],wn[k],Legendre(N,xn[k])))
        
    return xn,wn

#============================
def ChebyshevRoots(N):
    xn = np.linspace(-1,1,N)
    wn = np.linspace(-1,1,N)
    
    for n in range(0,N):
        ns=n+1
        xn[n]=math.cos((2*ns-1)*math.pi/2/N)
        wn[n]=math.pi/N
        
        #print("Finding root %d of %d, root=%f, weight=%f" %(n+1,N,xn[n],wn[n]))
        
    return xn,wn
           
#============================
def AssociatedLegendre(ell,m,x):
    if (abs(m)>ell):
        return 0.0

    if ell==0:
        return 1.0
    
    #====m=0,l=1
    Pn = x
    
    #====m=1,l=1
    Pnpos= -math.sqrt(1-x**2)

    #====m=-1,l=1
    Pnneg= -0.5*Pnpos
    
    if (ell==1):
        if (m==-1):
            return Pnneg
        if (m==0):
            return Pn
        if (m==1):
            return Pnpos

    Pmlp1 = 0
    if (ell==m):
        Pmlp1 = -(2.0*ell-1)*math.sqrt(1-x**2) * \
                AssociatedLegendre(ell-1,ell-1,x)
    else:
        Pmlp1 = (2.0*ell-1)*x*AssociatedLegendre(ell-1,m,x)
        Pmlp1 = Pmlp1 - (ell+m-1)*AssociatedLegendre(ell-2,m,x)
        Pmlp1 = Pmlp1 / (ell-m)
    
    return Pmlp1 

#============================
def fac(x):
    if (x==0):
        return 1
    if (x==1):
        return 1
    
    return fac(x-1)*x

#============================
def Min1powerm(m):
    if (m==0):
        return 1
    if ((m%2)==0):
        return 1
    else:
        return -1

#============================    
def Ylm(ell,m,varphi,theta):
    if (m<0):
        return Min1powerm(abs(m))*math.sqrt( \
               2.0 * \
               fac(ell-abs(m))/fac(ell+abs(m)) )* \
               AssociatedLegendre(ell,abs(m),math.cos(theta))* \
               math.sin(abs(m)*varphi)
    elif (m==0):
        return AssociatedLegendre(ell,abs(m),math.cos(theta))
    else:
        return Min1powerm(abs(m))*math.sqrt( \
               2.0 * \
               fac(ell-abs(m))/fac(ell+abs(m)) )* \
               AssociatedLegendre(ell,abs(m),math.cos(theta))* \
               math.cos(abs(m)*varphi)

#============================
class Quadrature:

    #===================
    def InitializeWithGL(self,Np):
        self.Na = 1
        self.Np = Np*2

        self.xp,self.wp = LegendreRoots(self.Np)
        self.wa = np.array([2.0*math.pi])

        self.theta = np.zeros(self.Np)
        for i in range(0, self.Np):
            self.theta[i] = math.acos(self.xp[i])

        self.varphi = np.array([0.0])

        self.omegas = []
        self.weights = []
        for i in range(0, self.Na):
            varphi = self.varphi[i]
            for j in range(0, self.Np):
                theta = self.theta[j]
                omega = np.array([math.sin(theta) * math.cos(varphi),
                                  math.sin(theta) * math.sin(varphi),
                                  math.cos(theta)])
                self.omegas.append(omega)
                self.weights.append(self.wp[j]*self.wa[i])

        
    #===================
    def InitializeWithGLC(self,Na,Np):
        self.Na = Na*4
        self.Np = Np*2
        self.xp,self.wp = LegendreRoots(self.Np)
        self.xa,self.wa = ChebyshevRoots(self.Na)
        
        self.wa=2*self.wa
        
        self.theta = np.zeros(self.Np)
        for i in range(0,self.Np):
            self.theta[i] = math.acos(self.xp[i])
        
        self.varphi = np.zeros(self.Na)
        for i in range(0,self.Na):
            self.varphi[i] = math.pi*(2*(i+1)-1)/self.Na

        self.omegas = []
        self.weights = []
        for i in range(0, self.Na):
            varphi = self.varphi[i]
            for j in range(0, self.Np):
                theta = self.theta[j]
                omega = np.array([math.sin(theta) * math.cos(varphi),
                                  math.sin(theta) * math.sin(varphi),
                                  math.cos(theta)])
                self.omegas.append(omega)
                self.weights.append(self.wp[j]*self.wa[i])
       
    #===================
    def InitializeWithGLL(self,Na,Np):
        self.Na = Na*4
        self.Np = Np*2
        self.xp,self.wp = LegendreRoots(self.Np)
        self.xa,self.wa = LegendreRoots(self.Na)
        
        self.wa=math.pi*self.wa
        
        self.theta = np.zeros(self.Np)
        for i in range(0,self.Np):
            self.theta[i] = math.acos(self.xp[i])
        
        self.varphi = np.zeros(self.Na)
        for i in range(0,self.Na):
            self.varphi[i] = math.pi*self.xa[i]+math.pi

        self.omegas = []
        self.weights = []
        for i in range(0, self.Na):
            varphi = self.varphi[i]
            for j in range(0, self.Np):
                theta = self.theta[j]
                omega = np.array([math.sin(theta) * math.cos(varphi),
                                  math.sin(theta) * math.sin(varphi),
                                  math.cos(theta)])
                self.omegas.append(omega)
                self.weights.append(self.wp[j]*self.wa[i])