'''

BASIC FUNCTIONS

Keep functions in alphabeticl order.

'''


# Required packages
import numpy as np
from matplotlib import path
from scipy import stats, linalg
import pickle
from scipy.io import loadmat

'''============================================================================
PICKLE 
============================================================================'''

def ipickle(path):
    fin=open(path,'rb')
    data=pickle.load(fin)
    fin.close()
    return data

def opickle(path,data):
    fout=open(path,'wb')
    pickle.dump(data,fout); 
    fout.close()

'''============================================================================
ADD LETTERS TO FIGURE PANELS
============================================================================'''

def axletters(ax,plt,rx,ry,*args):
    lab=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)']
    labCap=['A','B','C','D','E','F','G','H','I','J']
    c=0
    if len(ax.shape)==2:
        for i in range(ax.shape[0]):
            for j in range(ax.shape[1]):
                yl=ax[i,j].get_ylim()
                xl=ax[i,j].get_xlim()
                dyl=np.abs(yl[1]-yl[0])
                dxl=np.abs(xl[1]-xl[0])
                ax[i,j].text(xl[0]+rx*dxl,yl[0]+ry*dyl,lab[c],
                  {'color':plt.rcParams.get('axes.edgecolor'),
                   'fontsize':plt.rcParams.get('font.size')})
                try:
                    ax[i,j].text(xl[0]+(rx+args[1])*dxl,yl[0]+ry*dyl,args[0][c],
                      {'color':plt.rcParams.get('axes.edgecolor'),
                       'fontsize':plt.rcParams.get('font.size')})
                except:
                    pass
                    
                c=c+1
    elif len(ax.shape)==1:
        for i in range(ax.shape[0]):
            yl=ax[i].get_ylim()
            xl=ax[i].get_xlim()
            dyl=np.abs(yl[1]-yl[0])
            dxl=np.abs(xl[1]-xl[0])
            ax[i].text(xl[0]+rx*dxl,yl[0]+ry*dyl,lab[c],
              {'color':plt.rcParams.get('axes.edgecolor'),
               'fontsize':plt.rcParams.get('font.size')})
            try:
                ax[i].text(xl[0]+(rx+args[1])*dxl+args[1],yl[0]+ry*dyl,args[0][c],
                  {'color':plt.rcParams.get('axes.edgecolor'),
                   'fontsize':plt.rcParams.get('font.size')})
            except:
                pass
            c=c+1

'''============================================================================
BUNCH CONTENTS OF DICTIONARY
============================================================================'''

class Bunch(object):
    def __init__(self, adict):
        self.__dict__.update(adict)

'''============================================================================
CONVERT CM TO INCHES
============================================================================'''

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
    
'''============================================================================
COUNT BY CATEGORIES
============================================================================'''

def CountByCategories(z,*args):

    z=z.flatten()
    u=np.unique(z)
    N=np.zeros(u.size)
    for i in range(u.size):
        ind=np.where(z==u[i])[0]
        N[i]=ind.size

    try:
        if args[0]=='Percent':
            N=N/np.sum(N)*100
    except:
        pass

    N=np.round(N,2)
    return u,N

'''============================================================================
ISEMPTY
============================================================================'''

def isempty(value):
    try:
        value = float(value)
    except ValueError:
        pass
    return bool(value)

'''============================================================================
IMPORT MATLAB .MAT FILE
============================================================================'''

def ImportMat(pth,vnam):
    mat=loadmat(pth)
    d={n: mat[vnam][n][0, 0] for n in mat[vnam].dtype.names}
    return d

'''============================================================================
MINIMUM AND MAXIMUM
============================================================================'''

def minmax(ar):
    
    mm=[np.amin(ar),np.amax(ar)]
    
    return mm
    


'''============================================================================
MOVING AVERAGE
============================================================================'''

def movingave(y0,period,meth):
    
    if y0.ndim==1:
        # Convert to 2d array
        y=y0.reshape((y0.shape[0],1))
    else:
        y=y0
    
    m,n=y.shape   
    iy=np.arange(0,m,1)    
    ma=np.nan*np.ones((m,n))
    for i in range(m):
        if (meth=='Centre') | (meth=='centre') | (meth=='center') | (meth=='Center'):
            ind=np.where( (iy>=np.round(i-period/2)) & (iy<=np.round(i+period/2)) )[0]
        elif (meth=='Historical') | (meth=='historical'):
            ind=np.where( (iy>=np.round(i-period)) & (iy<=i) )[0]
        mu=np.nanmean(y[ind,:],axis=0)
        ma[i,:]=mu
    
    if y0.ndim==1:
        ma=ma.flatten()
        
    return ma
 
'''============================================================================
PARTIAL CORRELATION
============================================================================'''

def PartialCorrelation(C):
    C=np.asarray(C)
    p=C.shape[1]
    P_corr=np.zeros((p, p), dtype=np.float)
    for i in range(p):
        P_corr[i, i]=1
        for j in range(i+1, p):
            idx=np.ones(p, dtype=np.bool)
            idx[i]=False
            idx[j]=False
            beta_i=linalg.lstsq(C[:, idx], C[:, j])[0]
            beta_j=linalg.lstsq(C[:, idx], C[:, i])[0]

            res_j=C[:, j] - C[:, idx].dot( beta_i)
            res_i=C[:, i] - C[:, idx].dot(beta_j)

            corr=stats.pearsonr(res_i, res_j)[0]
            P_corr[i, j]=corr
            P_corr[j, i]=corr
    return P_corr
     
'''============================================================================
TIME VECTOR
============================================================================'''

def tvec(res,year1,year2):
    
    yr=np.arange(year1,year2,1)
    
    if res=='m':
        
        tv=np.zeros((yr.size*12,2),dtype=int)
        cnt=0
        for i in range(yr.size): 
            for j in range(12):
                tv[cnt,0]=yr[i]
                tv[cnt,1]=j+1
                cnt=cnt+1
            
        tv.shape    
                
    return tv





