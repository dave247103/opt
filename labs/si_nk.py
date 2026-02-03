import numpy as np, os
from io import StringIO
from plask import material

p=os.path.join(os.path.dirname(__file__),'Schinke.csv')
t=[s.strip() for s in open(p) if s.strip()]
i=next(j for j,s in enumerate(t) if s.lower()=='wl,k')
a=np.loadtxt(StringIO('\n'.join(t[1:i])),delimiter=',')
b=np.loadtxt(StringIO('\n'.join(t[i+1:])),delimiter=',')
wl=a[:,0]*1000
nr=a[:,1]
k=b[:,1]

@material.simple()
class Si_Shinke(material.Material):
    def Nr(self,lam,T=300.,n=0.): x=float(lam); return np.interp(x,wl,nr)+1j*np.interp(x,wl,k)
    def nr(self,lam,T=300.,n=0.): return float(np.interp(float(lam),wl,nr))
    def absp(self,lam,T=300.): x=float(lam); return 4*np.pi*np.interp(x,wl,k)*1e7/x