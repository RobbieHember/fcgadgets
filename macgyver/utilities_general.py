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
from scipy import stats
from subprocess import call
import netCDF4 as nc
import statsmodels.api as sm
import csv

#%% Regression stats

def GetRegStats(x,y):

    x=sm.tools.tools.add_constant(x)
    md=sm.OLS(y,x).fit()
    xhat=np.linspace(np.min(x[:,1]),np.max(x[:,1]),10)
    yhat=md.predict(np.c_[np.ones(xhat.size),xhat])

    rs={}
    rs['B']=md.params
    rs['pValues']=md.pvalues
    rs['R2']=md.rsquared
    rs['R2 adj']=md.rsquared_adj
    rs['AIC']=md.aic
    rs['P']=md.f_pvalue
    rs['xhat']=xhat
    rs['yhat']=yhat

    if rs['B'][0]<0:
        fl=' - '
    else:
        fl=' + '

    txt='y = ' + str(np.round(rs['B'][1],decimals=3)) + 'x' + fl + str(np.round(np.abs(rs['B'][0]),decimals=3)) + '\n' + 'R2 = ' + str(np.round(rs['R2'],decimals=2)) + '\n' + 'p < ' + str(np.round(rs['P'],decimals=2))
    rs['txt']=txt

    return rs,txt

#%% Import CSV file and add to dictionary

def ReadCSV(fin):
    d=pd.read_csv(fin).to_dict('list')
    for k in d.keys():
        d[k]=np.array(d[k])
    return d

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
    return

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
        df=pd.read_excel(args[0],sheet_name=args[1],skiprows=args[2])

    #df=df.where(pd.notnull(df),None)

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

def axletters(ax,plt,rx,ry,**kwargs):

    # Letter style
    Letter=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)']
    if 'LetterStyle' in kwargs.keys():
        if kwargs['LetterStyle']=='Caps':
            Letter=['A','B','C','D','E','F','G','H','I','J']
    # Font size
    fs=plt.rcParams.get('font.size')
    if 'FontSize' in kwargs.keys():
        fs=kwargs['FontSize']

    # Font weight
    fw='normal'
    if 'FontWeight' in kwargs.keys():
        if (kwargs['FontWeight']=='Bold') | (kwargs['FontWeight']=='bold'):
            fw='bold'

    # Font color
    fcl=plt.rcParams.get('axes.edgecolor')
    if 'FontColor' in kwargs.keys():
        print(kwargs['FontColor'])
        fcl=kwargs['FontColor']

    # Additional labels
    Label=['','','','','','','','','','','','','','']
    if 'Labels' in kwargs.keys():
        Label=kwargs['Labels']

    # Label spacer
    LabelSpacer=0.05
    if 'LabelSpacer' in kwargs.keys():
        LabelSpacer=kwargs['LabelSpacer']

    c=0
    if len(ax.shape)>1:
        for i in range(ax.shape[0]):
            for j in range(ax.shape[1]):

                yl=ax[i,j].get_ylim()
                xl=ax[i,j].get_xlim()

                dyl=np.abs(yl[1]-yl[0])
                dxl=np.abs(xl[1]-xl[0])

                # Letter
                x=xl[0]+rx*dxl
                y=yl[0]+ry*dyl
                ax[i,j].text(x,y,Letter[c],{'color':fcl,'fontsize':fs,'fontweight':fw})

                # Label
                if Label[c]!='':
                    x=xl[0]+(rx+LabelSpacer)*dxl
                    y=yl[0]+ry*dyl
                    txt=Label[c]
                    ax[i,j].text(x+LabelSpacer,y,txt,{'color':fcl,'fontsize':fs,'fontweight':'normal'})

                c=c+1

    elif len(ax.shape)==1:

        for i in range(ax.shape[0]):

            yl=ax[i].get_ylim()
            xl=ax[i].get_xlim()

            dyl=np.abs(yl[1]-yl[0])
            dxl=np.abs(xl[1]-xl[0])

            # Letter
            x=xl[0]+rx*dxl
            y=yl[0]+ry*dyl
            ax[i].text(x,y,Letter[c],{'color':fcl,'fontsize':fs,'fontweight':fw})

            # Label
            if Label[c]!='':
                x=xl[0]+(rx+LabelSpacer)*dxl
                y=yl[0]+ry*dyl
                txt=Label[c]
                ax[i].text(x,y,txt,{'color':fcl,'fontsize':fs,'fontweight':'normal'})

            c=c+1

    return

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

#%% INTERSECT ARRAYS (RETURNING ALL INDICES)
# *** doesn't work! ***

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

    yr=np.arange(year1,year2+1,1)

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

# N,mu,med,sig,se=gu.discres(x,y,bw,bin)

def discres(x,y,bw,bin):
    N=np.nan*np.ones(bin.size)
    mu=np.nan*np.ones(bin.size)
    med=np.nan*np.ones(bin.size)
    sig=np.nan*np.ones(bin.size)
    se=np.nan*np.ones(bin.size)
    for i in range(bin.size):
        ind=np.where(np.abs(x-bin[i])<=bw/2)[0]
        N[i]=ind.size
        if ind.size==0:
            continue
        mu[i]=np.nanmean(y[ind])
        med[i]=np.nanmedian(y[ind])
        sig[i]=np.nanstd(y[ind])
        se[i]=np.nanstd(y[ind])/np.sqrt(ind.size)
    return N,mu,med,sig,se

#%% Binned sum

def BinnedCount(bins,bw,x):

    y=np.zeros(bins.size)
    for i in range(bins.size):
        ind=np.where(np.abs(x-bins[i])<=bw/2)[0]
        y[i]=ind.size

    return y

#%% GRAPHICS PARAMETERS FOR JUPYTER NOTEBOOK

# def Import_GraphicsParameters(type):
#     if type=='ipynb':
#         fs1=6
#         fs2=9
#         params={'font.sans-serif':'Arial',
#                 'font.size':fs1,
#                 'figure.titlesize':fs2,
#                 'figure.dpi':150,
#                 'figure.constrained_layout.use':True,
#                 'axes.edgecolor':'black',
#                 'axes.labelsize':fs1,
#                 'axes.labelcolor':'black',
#                 'axes.titlesize':fs2,
#                 'axes.titlepad':2,
#                 'axes.linewidth':0.5,
#                 'lines.linewidth':0.75,
#                 'lines.markersize': 3.0,
#                 'text.color':'black',
#                 'xtick.color':'black',
#                 'xtick.labelsize':fs1,
#                 'xtick.major.width':0.5,
#                 'xtick.major.size':3,
#                 'xtick.direction':'in',
#                 'ytick.color':'black',
#                 'ytick.labelsize':fs1,
#                 'ytick.major.width':0.5,
#                 'ytick.major.size':3,
#                 'ytick.direction':'in',
#                 'legend.fontsize':fs1,
#                 'savefig.dpi':900,
#                 'savefig.transparent':True,
#                 'savefig.format':'png',
#                 'savefig.pad_inches':0.1,
#                 'savefig.bbox':'tight'}
#     elif type=='spyder_fs6':
#         fs1=6
#         fs2=6
#         params={'font.sans-serif':'Arial',
#                 'font.size':fs1,
#                 'figure.titlesize':fs2,
#                 'figure.dpi':150,
#                 'figure.constrained_layout.use':True,
#                 'axes.edgecolor':'black',
#                 'axes.labelsize':fs1,
#                 'axes.labelcolor':'black',
#                 'axes.titlesize':fs2,
#                 'axes.titlepad':2,
#                 'axes.linewidth':0.5,
#                 'lines.linewidth':0.75,
#                 'lines.markersize': 2,
#                 'text.color':'black',
#                 'xtick.color':'black',
#                 'xtick.labelsize':fs1,
#                 'xtick.major.width':0.5,
#                 'xtick.major.size':3,
#                 'xtick.direction':'in',
#                 'ytick.color':'black',
#                 'ytick.labelsize':fs1,
#                 'ytick.major.width':0.5,
#                 'ytick.major.size':3,
#                 'ytick.direction':'in',
#                 'legend.fontsize':fs1,
#                 'savefig.dpi':900,
#                 'savefig.transparent':True,
#                 'savefig.format':'png',
#                 'savefig.pad_inches':0.1,
#                 'savefig.bbox':'tight'}
#     elif type=='spyder_fs7':
#         fs1=7
#         fs2=7
#         params={'font.sans-serif':'Arial',
#                 'font.size':fs1,
#                 'figure.titlesize':fs2,
#                 'figure.dpi':150,
#                 'figure.constrained_layout.use':True,
#                 'axes.edgecolor':'black',
#                 'axes.labelsize':fs1,
#                 'axes.labelcolor':'black',
#                 'axes.titlesize':fs2,
#                 'axes.titlepad':2,
#                 'axes.linewidth':0.5,
#                 'lines.linewidth':0.75,
#                 'lines.markersize': 2,
#                 'text.color':'black',
#                 'xtick.color':'black',
#                 'xtick.labelsize':fs1,
#                 'xtick.major.width':0.5,
#                 'xtick.major.size':3,
#                 'xtick.direction':'in',
#                 'ytick.color':'black',
#                 'ytick.labelsize':fs1,
#                 'ytick.major.width':0.5,
#                 'ytick.major.size':3,
#                 'ytick.direction':'in',
#                 'legend.fontsize':fs1,
#                 'savefig.dpi':900,
#                 'savefig.transparent':False,
#                 'savefig.format':'png',
#                 'savefig.pad_inches':0.1,
#                 'savefig.bbox':'tight'}
#     elif type=='Presentation':
#         fs1=7
#         fs2=7
#         params={'font.sans-serif':'Arial',
#                 'font.size':fs1,
#                 'figure.titlesize':fs2,
#                 'figure.dpi':150,
#                 'figure.constrained_layout.use':True,
#                 'axes.edgecolor':'black',
#                 'axes.labelsize':fs1,
#                 'axes.labelcolor':'black',
#                 'axes.titlesize':fs2,
#                 'axes.titlepad':2,
#                 'axes.linewidth':0.5,
#                 'lines.linewidth':1.25,
#                 'lines.markersize': 4,
#                 'text.color':'black',
#                 'xtick.color':'black',
#                 'xtick.labelsize':fs1,
#                 'xtick.major.width':0.5,
#                 'xtick.major.size':3,
#                 'xtick.direction':'in',
#                 'ytick.color':'black',
#                 'ytick.labelsize':fs1,
#                 'ytick.major.width':0.5,
#                 'ytick.major.size':3,
#                 'ytick.direction':'in',
#                 'legend.fontsize':fs1,
#                 'savefig.dpi':900,
#                 'savefig.transparent':True,
#                 'savefig.format':'png',
#                 'savefig.pad_inches':0.1,
#                 'savefig.bbox':'tight'}
#     elif type=='Presentation Dark':

#         fs1=7; fs2=7; cla=[0.92,0.92,0.92];
#         params={'font.sans-serif':'Arial',
#                 'font.size':fs1,
#                 'figure.titlesize':fs2,
#                 'figure.dpi':150,
#                 'figure.constrained_layout.use':True,
#                 'axes.edgecolor':cla,
#                 'axes.labelsize':fs1,
#                 'axes.labelcolor':cla,
#                 'axes.titlesize':fs2,
#                 'axes.titlepad':2,
#                 'axes.linewidth':0.5,
#                 'lines.linewidth':5,
#                 'lines.markersize': 4,
#                 'text.color':cla,
#                 'xtick.color':cla,
#                 'xtick.labelsize':fs1,
#                 'xtick.major.width':0.5,
#                 'xtick.major.size':3,
#                 'xtick.direction':'in',
#                 'ytick.color':cla,
#                 'ytick.labelsize':fs1,
#                 'ytick.major.width':0.5,
#                 'ytick.major.size':3,
#                 'ytick.direction':'in',
#                 'legend.fontsize':fs1,
#                 'savefig.dpi':900,
#                 'savefig.transparent':True,
#                 'savefig.format':'png',
#                 'savefig.pad_inches':0.1,
#                 'savefig.bbox':'tight'}

#     return params

#%% Set graphics
def SetGraphics(type):

    gp={}

    if type=='Presentation Dark':

        gp['fs1']=7
        gp['fs2']=7
        gp['fs3']=6
        gp['cla']=[0.9,0.9,0.9]
        gp['clt']=[0.8,0.8,0.8]
        gp['cl1']=[0.27,0.49,0.79];
        gp['cl2']=[0.6,0.9,0];
        gp['cl3']=[0.7,0.4,0.95];
        gp['cl4']=[0.85,0,0]
        gp['lw1']=0.5;
        gp['lw2']=1.0;
        gp['ms']=4
        gp['Alpha1']=0.225;
        gp['Alpha2']=0.45;
        gp['tickl']=1.5;

        params={'font.sans-serif':'Arial',
                'font.size':gp['fs1'],
                'figure.titlesize':gp['fs2'],
                'figure.constrained_layout.use':True,
                'axes.edgecolor':gp['cla'],
                'axes.labelsize':gp['fs1'],
                'axes.labelcolor':gp['cla'],
                'axes.titlesize':gp['fs2'],
                'axes.titlepad':2,
                'axes.linewidth':gp['lw1'],
                'lines.linewidth':gp['lw2'],
                'lines.markersize':gp['ms'],
                'text.color':gp['cla'],
                'xtick.color':gp['cla'],
                'xtick.labelsize':gp['fs1'],
                'xtick.major.width':gp['lw1'],
                'xtick.major.size':3,
                'xtick.direction':'in',
                'ytick.color':gp['cla'],
                'ytick.labelsize':gp['fs1'],
                'ytick.major.width':0.5,
                'ytick.major.size':3,
                'ytick.direction':'in',
                'legend.fontsize':gp['fs1'],
                'savefig.dpi':900,
                'savefig.transparent':True,
                'savefig.format':'png',
                'savefig.pad_inches':0.1,
                'savefig.bbox':'tight'}

    elif type=='Presentation Light':

        gp['fs1']=9
        gp['fs2']=7
        gp['fs3']=6
        gp['cla']=[0.5,0.5,0.5]
        gp['clt']=[0.5,0.5,0.5]
        gp['cl1']=[0.27,0.49,0.79];
        gp['cl2']=[0.6,0.9,0];
        gp['cl3']=[0.7,0.4,0.95];
        gp['cl4']=[0.85,0,0]
        gp['lw1']=0.5;
        gp['lw2']=1.0;
        gp['ms']=4
        gp['Alpha1']=0.225;
        gp['Alpha2']=0.45;
        gp['tickl']=1.5;

        params={'font.sans-serif':'Arial',
                'font.size':gp['fs1'],
                'figure.titlesize':gp['fs2'],
                'figure.constrained_layout.use':True,
                'axes.edgecolor':gp['cla'],
                'axes.labelsize':gp['fs1'],
                'axes.labelcolor':gp['cla'],
                'axes.titlesize':gp['fs2'],
                'axes.titlepad':2,
                'axes.linewidth':gp['lw1'],
                'lines.linewidth':gp['lw2'],
                'lines.markersize':gp['ms'],
                'text.color':gp['cla'],
                'xtick.color':gp['cla'],
                'xtick.labelsize':gp['fs1'],
                'xtick.major.width':gp['lw1'],
                'xtick.major.size':3,
                'xtick.direction':'in',
                'ytick.color':gp['cla'],
                'ytick.labelsize':gp['fs1'],
                'ytick.major.width':0.5,
                'ytick.major.size':3,
                'ytick.direction':'in',
                'legend.fontsize':gp['fs1'],
                'savefig.dpi':900,
                'savefig.transparent':True,
                'savefig.format':'png',
                'savefig.pad_inches':0.1,
                'savefig.bbox':'tight'}

    elif type=='Manuscript':

        gp['fs1']=7
        gp['fs2']=6
        gp['fs3']=5
        gp['cla']=[0,0,0]
        gp['clt']=[0,0,0]
        gp['cl1']=[0.27,0.49,0.79];
        gp['cl2']=[0.6,0.9,0];
        gp['cl3']=[0.7,0.4,0.95];
        gp['cl4']=[0.85,0,0]
        gp['lw1']=0.5;
        gp['lw2']=1.0;
        gp['ms']=4
        gp['Alpha1']=0.225;
        gp['Alpha2']=0.45;
        gp['tickl']=1.5;

        params={'font.sans-serif':'Arial',
                'font.size':gp['fs1'],
                'figure.titlesize':gp['fs2'],
                'figure.constrained_layout.use':True,
                'axes.edgecolor':gp['cla'],
                'axes.labelsize':gp['fs1'],
                'axes.labelcolor':gp['cla'],
                'axes.titlesize':gp['fs2'],
                'axes.titlepad':2,
                'axes.linewidth':gp['lw1'],
                'lines.linewidth':gp['lw2'],
                'lines.markersize':gp['ms'],
                'text.color':gp['cla'],
                'xtick.color':gp['cla'],
                'xtick.labelsize':gp['fs1'],
                'xtick.major.width':gp['lw1'],
                'xtick.major.size':3,
                'xtick.direction':'in',
                'ytick.color':gp['cla'],
                'ytick.labelsize':gp['fs1'],
                'ytick.major.width':0.5,
                'ytick.major.size':3,
                'ytick.direction':'in',
                'legend.fontsize':gp['fs1'],
                'savefig.dpi':900,
                'savefig.transparent':False,
                'savefig.format':'png',
                'savefig.pad_inches':0.1,
                'savefig.bbox':'tight'}

    plt.rcParams.update(params)

    return gp

#%% Z-score a variable

def zscore(x):

    if x.ndim==1:
        mu=np.nanmean(x)
        sig=np.nanstd(x)
        y=(x-mu)/sig
    else:
        y=x.copy()
        mu=np.zeros(x.shape[1])
        sig=np.zeros(x.shape[1])
        for j in range(x.shape[1]):
            mu[j]=np.nanmean(x[:,j])
            sig[j]=np.nanstd(x[:,j])
            y[:,j]=(x[:,j]-mu[j])/sig[j]

    return y,mu,sig

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
        c=np.arange(0,x.size,ivl,dtype=int)
        y=np.zeros(c.size)
        for i in range(c.size):
            y[i]=np.mean(x[c[i]:c[i]+ivl])
    else:
        y='Oops!'
        print('y=Oops, needs revision')
    return y

#%% Get confidnece intervals for the difference in two variables with confidence intervals

def GetCIsFromDifference(b_lo,b_hi,p_lo,p_hi):
    tmp=np.array([p_lo-b_lo,p_lo-b_hi,p_hi-b_lo,p_hi-b_hi]).T
    try:
        lo=np.min(tmp,axis=1)
        hi=np.max(tmp,axis=1)
    except:
        lo=np.min(tmp,axis=0)
        hi=np.max(tmp,axis=0)

    return lo,hi

#%% Get clamp

def Clamp(x,xmin,xmax):

    if x.size==1:
        y=np.maximum(np.minimum(xmax,x),xmin)
    else:
        y=np.max(np.min(xmax,x),xmin)

    return y

#%% Density profile

def ksdensity(y):
    kde=stats.gaussian_kde(y)
    x=np.linspace(y.min(),y.max(),100)
    p=kde(x)
    return p

#%% BUNCH CONTENTS OF DICTIONARY

class Bunch(object):
    def __init__(self, adict):
        self.__dict__.update(adict)
