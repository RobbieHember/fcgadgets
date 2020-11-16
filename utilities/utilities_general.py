'''

GENERAL FUNCTIONS

'''


#%% Import packages

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib import path
from scipy import stats, linalg
import pickle
from scipy.io import loadmat
from subprocess import call
import netCDF4 as nc

#%% PICKLE INPUT AND OUTPUT

def ipickle(path):
    fin=open(path,'rb')
    data=pickle.load(fin)
    fin.close()
    return data

def opickle(path,data):
    fout=open(path,'wb')
    pickle.dump(data,fout); 
    fout.close()

#%% PRINT FIGURE

def PrintFig(fin,type,dpi):    
    if type=='emf':
        plt.savefig(fin+'.svg',format='svg')
        try:
            path_to_inkscape=r'C:\Program Files\Inkscape\inkscape.exe'  
            call([path_to_inkscape,'--file',fin+'.svg','--export-emf',fin+'.emf'])            
        except:
            path_to_inkscape=r'C:\Program Files\Inkscape\bin\inkscape.exe'    
            call([path_to_inkscape,'--file',fin+'.svg','--export-emf',fin+'.emf'])
        os.remove(fin+'.svg')
    elif type=='png':
        plt.savefig(fin+'.png',format='png',dpi=dpi) 
    elif (type=='all') | (type=='All'):
        path_to_inkscape='C:\Program Files\Inkscape\inkscape.exe'    
        plt.savefig(fin+'.svg',format='svg')    
        call([path_to_inkscape,'--file',fin+'.svg','--export-emf',fin+'.emf'])
        os.remove(fin+'.svg')
        plt.savefig(fin+'.png',format='png',dpi=dpi) 
        
#%% IMPORT EXCEL SPREADSHEET AND CONVERT TO DICTIONARY

def ReadExcel(*args):
    
    if len(args)==1:
        df=pd.read_excel(args[0])
    elif len(args)==2:
        df=pd.read_excel(args[0],sheet_name=args[1])
    elif len(args)==3:
        df=pd.read_excel(args[0],skiprows=args[2])
    
    d=df.to_dict('list')
    for k in d.keys():
        d[k]=np.array(d[k])
    
    return d

#%% CONVERT DATAFRAME TO DICTIONARY

def DataFrameToDict(df):
    d=df.to_dict('list')
    for k in d.keys():
        d[k]=np.array(d[k])
    return d

#%% CONVERT DATAFRAME TO MATLAB-LIKE DATA STRUCTURE

def DataFrameToDataStruct(df):
    
    # Convert dataframe to dictionary
    d=DataFrameToDict(df)   
    
    # Make sure there are no spaces in the keys
    d={x.replace(' ','_'): v  
    for x, v in d.items()} 

    # Make sure there are no "." 
    d={x.replace('.','_'): v  
    for x,v in d.items()} 
    
    # Bunch the dictionary
    d=Bunch(d)
    
    return d

#%% ADD LETTERS TO FIGURE PANELS

def axletters(ax,plt,rx,ry,*args):
    lab=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)']
    labCap=['A','B','C','D','E','F','G','H','I','J']
    c=0
    if len(ax.shape)>1:
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

#%% BUNCH CONTENTS OF DICTIONARY

class Bunch(object):
    def __init__(self, adict):
        self.__dict__.update(adict)

#%% CONVERT CM TO INCHES

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)
    
#%% COUNT BY CATEGORIES

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

#%% IMPORT MATLAB .MAT FILE AND ADD TO DICTIONARY

def ImportMat(pth,vnam):
    mat=loadmat(pth)
    d={n: mat[vnam][n][0, 0] for n in mat[vnam].dtype.names}
    return d

#%% INTERSECT ARRAYS

def intersect(a,b):
    a1,ia=np.unique(a,return_index=True)
    b1,ib=np.unique(b,return_index=True)
    aux=np.concatenate((a1, b1))
    aux.sort()
    c=aux[:-1][aux[1:] == aux[:-1]]
    return c,ia[np.isin(a1,c)],ib[np.isin(b1,c)]

#%% ISEMPTY

def isempty(value):
    try:
        value=float(value)
    except ValueError:
        pass
    return bool(value)

#%% MINIMUM AND MAXIMUM

def minmax(ar):
    
    mm=[np.amin(ar),np.amax(ar)]
    
    return mm

#%% MOVING AVERAGE

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
 
#%% PARTIAL CORRELATION

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
     
#%% TIME VECTOR

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

#%% DIESCRETE RESPONSE (BINNED RESPONSE)

def discres(x,y,bw,bin):
    N=np.nan*np.ones(bin.size)
    mu=np.nan*np.ones(bin.size)
    sig=np.nan*np.ones(bin.size)
    se=np.nan*np.ones(bin.size)
    for i in range(bin.size):
        ind=np.where(np.abs(x-bin[i])<=bw/2)[0]
        N[i]=ind.size
        if ind.size==0:
            continue
        mu[i]=np.nanmean(y[ind])
        sig[i]=np.nanstd(y[ind])
        se[i]=np.nanstd(y[ind])/np.sqrt(ind.size)    
    return N,mu,sig,se

#%% GRAPHICS PARAMETERS FOR JUPYTER NOTEBOOK

def Import_GraphicsParameters(type):
    if type=='ipynb':
        fs1=6
        fs2=9
        params={'font.sans-serif':'Arial',
                'font.size':fs1,
                'figure.titlesize':fs2,
                'figure.dpi':150,
                'figure.constrained_layout.use':True,
                'axes.edgecolor':'black',
                'axes.labelsize':fs1,
                'axes.labelcolor':'black',
                'axes.titlesize':fs2,
                'axes.titlepad':2,
                'axes.linewidth':0.5,        
                'lines.linewidth':0.75,
                'lines.markersize': 3.0,
                'text.color':'black',
                'xtick.color':'black',        
                'xtick.labelsize':fs1,
                'xtick.major.width':0.5,
                'xtick.major.size':3,
                'xtick.direction':'in',
                'ytick.color':'black',
                'ytick.labelsize':fs1,
                'ytick.major.width':0.5,
                'ytick.major.size':3,
                'ytick.direction':'in',
                'legend.fontsize':fs1,        
                'savefig.dpi':900,
                'savefig.transparent':True,
                'savefig.format':'png',
                'savefig.pad_inches':0.1,
                'savefig.bbox':'tight'}
    elif type=='spyder1':
        fs1=6
        fs2=7
        params={'font.sans-serif':'Arial',
                'font.size':fs1,
                'figure.titlesize':fs2,
                'figure.dpi':150,
                'figure.constrained_layout.use':True,
                'axes.edgecolor':'black',
                'axes.labelsize':fs1,
                'axes.labelcolor':'black',
                'axes.titlesize':fs2,
                'axes.titlepad':2,
                'axes.linewidth':0.5,        
                'lines.linewidth':0.75,
                'lines.markersize': 2,
                'text.color':'black',
                'xtick.color':'black',        
                'xtick.labelsize':fs1,
                'xtick.major.width':0.5,
                'xtick.major.size':3,
                'xtick.direction':'in',
                'ytick.color':'black',
                'ytick.labelsize':fs1,
                'ytick.major.width':0.5,
                'ytick.major.size':3,
                'ytick.direction':'in',
                'legend.fontsize':fs1,        
                'savefig.dpi':900,
                'savefig.transparent':True,
                'savefig.format':'png',
                'savefig.pad_inches':0.1,
                'savefig.bbox':'tight'}
    elif type=='spyder_fs7':
        fs1=7
        fs2=7
        params={'font.sans-serif':'Arial',
                'font.size':fs1,
                'figure.titlesize':fs2,
                'figure.dpi':150,
                'figure.constrained_layout.use':True,
                'axes.edgecolor':'black',
                'axes.labelsize':fs1,
                'axes.labelcolor':'black',
                'axes.titlesize':fs2,
                'axes.titlepad':2,
                'axes.linewidth':0.5,        
                'lines.linewidth':0.75,
                'lines.markersize': 2,
                'text.color':'black',
                'xtick.color':'black',        
                'xtick.labelsize':fs1,
                'xtick.major.width':0.5,
                'xtick.major.size':3,
                'xtick.direction':'in',
                'ytick.color':'black',
                'ytick.labelsize':fs1,
                'ytick.major.width':0.5,
                'ytick.major.size':3,
                'ytick.direction':'in',
                'legend.fontsize':fs1,        
                'savefig.dpi':900,
                'savefig.transparent':True,
                'savefig.format':'png',
                'savefig.pad_inches':0.1,
                'savefig.bbox':'tight'}
        
    return params

#%% Z-score a variable

def zscore(x):
    y=(x-np.nanmean(x))/np.nanstd(x)
    return y

#%% Read NetCDF4

def ReadNC(fin):
    ds=nc.Dataset(fin)
    d={}
    for k in ds.variables.keys():
        d[k]=np.array(ds.variables[k][:])    
    return d

#%% Sum over interval
    
def BlockSum(x,ivl):    
    if x.ndim==1:
        m=x.size
        x1=np.append(x,np.nan*np.ones(ivl))
        y=x[0::ivl]
        for i in range(1,ivl):
            y=y+x1[i::ivl]
    return y

#%% Sum over interval
    
def BlockMean(x,ivl):    
    if x.ndim==1:
        m=x.size
        y=np.zeros(np.ceil(x.size/ivl).astype(int))
        for i in range(ivl):
            to_add=x[i::ivl]
            y[0:to_add.size]=y[0:to_add.size]+to_add
        y=y/ivl
    else:
        y='Oops!'
        print('y=Oops, needs revision')        
    return y