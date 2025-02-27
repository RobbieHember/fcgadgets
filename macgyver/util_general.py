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
import calendar
import netCDF4 as nc
import statsmodels.api as sm
import csv

#%% Indices to repeated values in array
# idx=gu.IndicesFromUniqueArrayValues(x)
def IndicesFromUniqueArrayValues(x):
	if type(x[0])=='str':
		records_array=x.astype(str)
	else:
		records_array=x
	idx_sort=np.argsort(records_array)
	sorted_records_array=records_array[idx_sort]
	vals,idx_start,count=np.unique(sorted_records_array, return_counts=True, return_index=True)
	res=np.split(idx_sort,idx_start[1:])
	y=dict(zip(vals,res))
	return y

#%% Regression stats
# Eg: rs,txt=gu.GetRegStats(x,y)

# Scatterplot template:
flg=0
if flg==1:
	plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(7.8,7))
	ax.plot([-10000,10000],[-10000,10000],'k-',lw=2,color=[0.8,0.8,0.8])
	ax.plot(x,y,'o',ms=5,mec='w',mfc='k',mew=0.5)
	rs,txt=gu.GetRegStats(x,y)
	ax.plot(rs['xhat Line'],rs['yhat Line'],'r-')
	ax.text(1,1,rs['txt'],fontsize=7,ha='left')
	ax.text(20,20,'1:1',fontsize=7,ha='center')
	ax.set(xlabel='X lable',ylabel='Y label',xticks=np.arange(-30,30,5),yticks=np.arange(-30,30,5),xlim=[-30,30],ylim=[-30,30])
	ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
	plt.tight_layout()
	gu.PrintFig(meta['Graphics']['Print Figure Path'] + '\\Name','png',900)

def GetRegStats(x,y,*args):

	# Remove nans
	ikp=np.where(np.isnan(x+y)==False)[0]
	x=x[ikp]
	y=y[ikp]

	# Initialize data structure
	rs={}

	#--------------------------------------------------------------------------
	# Simple linear regression based on ordinary least squares
	#--------------------------------------------------------------------------
	
	x_WithOnes=sm.tools.tools.add_constant(x)
	if 'No Intercept' in args:
		md=sm.OLS(y,x).fit()
		yhat=md.predict(x)
	else:
		md=sm.OLS(y,x_WithOnes).fit()
		yhat=md.predict(x_WithOnes)
	
	rs['N']=x.size
	rs['DF']=x.size-md.params.size
	rs['B']=md.params
	rs['pValues']=md.pvalues
	rs['r']=np.sqrt(md.rsquared)
	rs['R2']=md.rsquared
	rs['R2 Adjusted']=md.rsquared_adj
	rs['AIC']=md.aic
	rs['P']=md.f_pvalue
	
	# Standard error of the estimate
	rs['Standard Error of the Estimate']=np.sqrt(np.sum((y-yhat)**2)/rs['N'])	
	
	# t-statistic for sample size N with degrees of freedom DF
	# (alpha = 0.05)
	t_two_tail=np.array([12.706,4.303,3.182,2.776,2.571,2.447,2.365,2.306,2.262,2.228,2.201,2.179,2.16,2.145,2.132,2.12,2.11,2.101,2.093,2.086,2.08,2.074,2.069,2.064,2.06,2.056,2.052,2.048,2.045,1.96])
	if rs['DF']<30:
	  t=t_two_tail[rs['DF']]
	else:
	  t=t_two_tail[-1]
	
	# 95% cofidence interval of the slope 
	rs['Slope Error']=1.96*np.sqrt(np.sum((y-yhat)**2)/rs['DF'])/np.sqrt(np.sum((x-np.mean(x))**2))
	
	# 95% cofidence interval of the y-intercept 
	rs['Intercept Error']=t*rs['Standard Error of the Estimate']*np.sqrt((1/rs['N'])+(np.mean(x)**2/(np.sum(x**2)-(np.sum(x)**2/rs['N']))))
	
	rs['xhat Line']=np.linspace(np.min(x),np.max(x),10)
	if 'No Intercept' in args:
		rs['yhat Line']=md.predict(rs['xhat Line'])
	else:
		rs['yhat Line']=md.predict(np.c_[np.ones(rs['xhat Line'].size),rs['xhat Line']])

	#--------------------------------------------------------------------------
	# Error between dependent and independent variables
	#--------------------------------------------------------------------------
	
	# Mean error (plus or minus S.E. of mean)
	me=np.mean(y)-np.mean(x)
	se_mean=2*(np.std(x-y)/np.sqrt(rs['N']))
	rs['Mean Error']=[me,se_mean]
	rs['Mean Error %']=me/np.mean(x)*100
	
	# Sum of squares
	rs['Sum of Squared Error']=np.sum((y-x)**2)
	
	# Root mean squared error
	rs['RMSE']=np.sqrt(rs['Sum of Squared Error']/rs['N'])
	rs['RMSE %']=rs['RMSE']/np.mean(x)*100
	
	# Mean absolute error
	rs['MAE']=np.mean(np.abs(y-x))
	rs['MAE %']=rs['MAE']/np.mean(x)*100
	
	# Mean squared error
	rs['MSE Total']=np.sum((x-y)**2)/rs['N']
	
	# Systematic and unsystematic squared errors (Willmott 1981)
	# The sum will equal the mean squared error
	rs['MSE Systematic']=(1/rs['N'])*np.sum((yhat-x)**2)	
	rs['MSE Unsystematic']=(1/rs['N'])*np.sum((y-yhat)**2)
	
	# Model efficiency coefficient (Nash and Sutcliffe 1970)
	rs['MEC']=1-np.sum((x-y)**2)/np.sum((x-np.mean(x))**2)

	#--------------------------------------------------------------------------
	# Text string
	#--------------------------------------------------------------------------
	
	if rs['B'][0]<0:
		fl=' - '
	else:
		fl=' + '

	if rs['P']<0.01:
		psmb='<'
	else:
		psmb='='

	if 'No Intercept' in args:
		txt='y = ' + str(np.round(rs['B'][0],decimals=3)) + 'x' + '\nR$^2$ = ' + str(np.round(rs['R2'],decimals=2)) + '\np ' + psmb + ' ' + str(np.round(rs['P'],decimals=2)) + '\nRMSE = ' + str(np.round(rs['RMSE'],decimals=2))
	else:
		txt='y = ' + str(np.round(rs['B'][1],decimals=3)) + 'x' + fl + str(np.round(np.abs(rs['B'][0]),decimals=3)) + '\nR$^2$ = ' + str(np.round(rs['R2'],decimals=2)) + '\np ' + psmb + ' ' + str(np.round(rs['P'],decimals=2)) + '\nRMSE = ' + str(np.round(rs['RMSE'],decimals=2))
	rs['txt']=txt

	return rs,txt

#%%
# rs,txt=gu.PolynomialFit(x,y,2)
def PolynomialFit(x,y,ord):

	# Remove nans
	ikp=np.where(np.isnan(x+y)==False)[0]
	x=x[ikp]
	y=y[ikp]

	x0=x.copy()

	# Initialize data structure
	rs={}

	if ord==2:
		x=np.column_stack((x,x**2))
	else:
		x=np.column_stack((x,x**2))

	x=sm.tools.tools.add_constant(x)
	md=sm.OLS(y,x).fit()
	#yhat=md.predict(x)
	
	rs['N']=x.shape[0]
	rs['DF']=x.size-md.params.size
	rs['B']=md.params
	rs['pValues']=md.pvalues
	rs['r']=np.sqrt(md.rsquared)
	rs['R2']=md.rsquared
	rs['R2 Adjusted']=md.rsquared_adj
	rs['AIC']=md.aic
	rs['P']=md.f_pvalue

	rs['xhat Line']=np.linspace(np.min(x0),np.max(x0),10)
	rs['yhat Line']=md.predict(np.c_[np.ones(rs['xhat Line'].size),rs['xhat Line'],rs['xhat Line']**2])

	#--------------------------------------------------------------------------
	# Text string
	#--------------------------------------------------------------------------
	
	if rs['B'][0]<0:
		fl=' - '
	else:
		fl=' + '

	if rs['P']<0.01:
		psmb='<'
	else:
		psmb='='


	txt=''
	rs['txt']=txt

	return rs,txt

#%% Import CSV file and add to dictionary
def ReadCSV(fin,skip_rows):
	if skip_rows==0:
		d=pd.read_csv(fin).to_dict('list')
	else:
		# *** USE THIS TO SKIP ROWS: ***
		d=pd.read_csv(fin,skiprows=[skip_rows],encoding='unicode_escape').to_dict('list')

	for k in d.keys():
		d[k]=np.array(d[k])

	return d

#%% PICKLE INPUT AND OUTPUT
def ipickle(path):
	data=pd.read_pickle(path)
	# Below started to give errors when opening old files
	#fin=open(path,'rb')
	#data=pickle.load(fin)
	#fin.close()
	return data

def opickle(path,data):
	fout=open(path,'wb')
	pickle.dump(data,fout);
	fout.close()
	return

#%%
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
		# *** this method of printing emf does not work, do it manually in Inkscape ***
		#path_to_inkscape=r'C:\Program Files\Inkscape\bin\inkscape.exe'
		plt.savefig(fin+'.svg',format='svg')
		#call([path_to_inkscape,'--file',fin+'.svg','--export-emf',fin+'.emf'])
		#os.remove(fin+'.svg')
		plt.savefig(fin+'.png',format='png',dpi=dpi)
	return

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

#%%
def ReadExcelToLine(pthin,n):
	df=pd.read_excel(pthin,nrows=n)
	d=df.to_dict('list')
	for k in d.keys():
		d[k]=np.array(d[k])
	return d

#%% Save dictionary to excel spreadsheet
# Example: gu.PrintDict(d,pth,SheetName='Sheet1')
def PrintDict(d,pth,**kwargs):

	if 'SheetName' in kwargs.keys():
		sn=kwargs['SheetName']
	else:
		sn='Sheet1'

	# Convert to dataframe
	df=pd.DataFrame.from_dict(d)

	# Export
	writer=pd.ExcelWriter(pth,engine="xlsxwriter")
	df.style.set_properties(**{'text-align':'left'}).to_excel(writer,sheet_name=sn,startrow=0,header=True,index=False)

	workbook=writer.book
	worksheet=writer.sheets[sn]
	
	header_format=workbook.add_format({'bold':True,'text_wrap':False,'valign':'top','align':'left','border':0}) # ,'fg_color':'#D7E4BC'
	
	for col_num,value in enumerate(df.columns.values):
		worksheet.write(0,col_num,value,header_format)

	for column in df:
		column_length=max(df[column].astype(str).map(len).max(),len(column))
		col_idx=df.columns.get_loc(column)
		writer.sheets[sn].set_column(col_idx,col_idx,column_length)

	writer.close()
	return

#%% Save dataframe to excel spreadsheet
# Example: gu.PrintDict(d,pth,SheetName='Sheet1')
def PrintDF(df,pth,**kwargs):

	if 'SheetName' in kwargs.keys():
		sn=kwargs['SheetName']
	else:
		sn='Sheet1'

	# Export
	writer=pd.ExcelWriter(pth,engine="xlsxwriter")
	df.style.set_properties(**{'text-align':'left'}).to_excel(writer,sheet_name=sn,startrow=0,header=True,index=False)

	workbook=writer.book
	worksheet=writer.sheets[sn]
	
	header_format=workbook.add_format({'bold':True,'text_wrap':False,'valign':'top','align':'left','border':0}) # ,'fg_color':'#D7E4BC'
	
	for col_num,value in enumerate(df.columns.values):
		worksheet.write(0,col_num,value,header_format)

	for column in df:
		column_length=max(df[column].astype(str).map(len).max(),len(column))
		col_idx=df.columns.get_loc(column)
		writer.sheets[sn].set_column(col_idx,col_idx,column_length)

	writer.close()
	return

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
	Letter=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)']
	if 'LetterStyle' in kwargs.keys():
		if kwargs['LetterStyle']=='Caps':
			Letter=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
		elif kwargs['LetterStyle']=='NoPar':
			Letter=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p']
		elif kwargs['LetterStyle']=='Default':
			Letter=['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)']
		else:
			pass

	# Font size
	fs=plt.rcParams.get('font.size')
	if 'FontSize' in kwargs.keys():
		fs=kwargs['FontSize'];

	# Font weight
	fw='normal'
	if 'FontWeight' in kwargs.keys():
		if (kwargs['FontWeight']=='Bold') | (kwargs['FontWeight']=='bold'):
			fw='bold'

	# Font color
	fcl=plt.rcParams.get('axes.edgecolor');
	if 'FontColor' in kwargs.keys():
		#print(kwargs['FontColor'])
		fcl=kwargs['FontColor'];

	# Skip axis
	Skip_Flag=np.zeros(ax.shape)
	if 'Skip' in kwargs.keys():
		Skip_Flag[kwargs['Skip']]=1

	# Additional labels
	Label=['','','','','','','','','','','','','','']
	if 'Labels' in kwargs.keys():
		Label=kwargs['Labels']

	# Label spacer
	LabelSpacer=0.05
	if 'LabelSpacer' in kwargs.keys():
		LabelSpacer=kwargs['LabelSpacer']

	# Start from non-A letter
	if 'StartLetterIndex' in kwargs.keys():
		c=kwargs['StartLetterIndex']
	else:
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

		for i in range(ax.size):
			
			if Skip_Flag[i]==1:
				continue

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

#%% Count By Category
# d=gu.CountByCategories(z,'Percent')
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

	d={}
	for i in range(u.size):
		d[u[i]]=N[i]

	return d

#%% Stats By Category
# sts=gu.StatsByCategories(vCat,vStat)
def StatsByCategories(vCat,vStat):
	sts={}
	sts['u']=np.unique(vCat)
	sts['N']=np.zeros(sts['u'].size)
	sts['mu']=np.zeros(sts['u'].size)
	sts['sum']=np.zeros(sts['u'].size)
	for i in range(sts['u'].size):
		ind=np.where(vCat==sts['u'][i])[0]
		sts['mu'][i]=np.nanmean(vStat[ind])
		sts['sum'][i]=np.nansum(vStat[ind])
		sts['N'][i]=ind.size
	return sts

#%% Stats By Category (Dictionary)
# sts=gu.DictStatsByCategories(vCat,lutCat,vStat)
def DictStatsByCategories(vCat,lutCat,vStat):
	sts={}
	sts['N']={}
	sts['mu']={}
	sts['sum']={}
	for k in lutCat.keys():
		ind=np.where(vCat==lutCat[k])[0]
		sts['mu'][k]=np.round(np.nanmean(vStat[ind]),decimals=2)
		sts['sum'][k]=np.round(np.nansum(vStat[ind]),decimals=2)
		sts['N'][k]=ind.size
	return sts

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

#%% Correlation coefficient (mimmicing MATLAB)
def corrcoef(dframe):
	fmatrix = dframe.values
	rows, cols = fmatrix.shape
	r = np.ones((cols, cols), dtype=float)
	p = np.ones((cols, cols), dtype=float)
	for i in range(cols):
		for j in range(cols):
			if i == j:
				r_, p_ = 1., 1.
			else:
				r_, p_ = stats.pearsonr(fmatrix[:,i], fmatrix[:,j])
			r[j][i] = r_
			p[j][i] = p_
	return r, p

#%%
def PartialCorrelation(C):
	C=np.asarray(C)
	p=C.shape[1]
	P_corr=np.zeros((p, p), dtype='float')
	for i in range(p):
		P_corr[i, i]=1
		for j in range(i+1, p):
			idx=np.ones(p, dtype='bool')
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

#%% Time vector
# tv=gu.tvec('m',1945,2100)
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
	elif res=='d':
		tv=np.zeros((1000000,4),dtype='int16')
		cnt=0
		for iY in range(yr.size):
			doy=1
			for iM in range(12):
				n=calendar.monthrange(yr[iY],iM+1)[1]
				for iD in range(n):
					tv[cnt,:]=np.array([yr[iY],iM+1,iD+1,doy])
					doy=doy+1
					cnt=cnt+1
		tv=tv[0:cnt,:]
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
		mu[i]=np.nanmean(y[ind]);
		med[i]=np.nanmedian(y[ind]);
		sig[i]=np.nanstd(y[ind]);
		se[i]=np.nanstd(y[ind])/np.sqrt(ind.size);
	return N,mu,med,sig,se

#%% Binned count
def BinnedCount(bins,bw,x):
	y=np.zeros(bins.size)
	for i in range(bins.size):
		ind=np.where(np.abs(x-bins[i])<=bw/2)[0]
		y[i]=ind.size
	return y

#%% GRAPHICS PARAMETERS FOR JUPYTER NOTEBOOK

# def Import_GraphicsParameters(type):
#	 if type=='ipynb':
#		 fs1=6
#		 fs2=9
#		 params={'font.sans-serif':'Arial',
#				 'font.size':fs1,
#				 'figure.titlesize':fs2,
#				 'figure.dpi':150,
#				 'figure.constrained_layout.use':True,
#				 'axes.edgecolor':'black',
#				 'axes.labelsize':fs1,
#				 'axes.labelcolor':'black',
#				 'axes.titlesize':fs2,
#				 'axes.titlepad':2,
#				 'axes.linewidth':0.5,
#				 'lines.linewidth':0.75,
#				 'lines.markersize': 3.0,
#				 'text.color':'black',
#				 'xtick.color':'black',
#				 'xtick.labelsize':fs1,
#				 'xtick.major.width':0.5,
#				 'xtick.major.size':3,
#				 'xtick.direction':'in',
#				 'ytick.color':'black',
#				 'ytick.labelsize':fs1,
#				 'ytick.major.width':0.5,
#				 'ytick.major.size':3,
#				 'ytick.direction':'in',
#				 'legend.fontsize':fs1,
#				 'savefig.dpi':900,
#				 'savefig.transparent':True,
#				 'savefig.format':'png',
#				 'savefig.pad_inches':0.1,
#				 'savefig.bbox':'tight'}
#	 elif type=='spyder_fs6':
#		 fs1=6
#		 fs2=6
#		 params={'font.sans-serif':'Arial',
#				 'font.size':fs1,
#				 'figure.titlesize':fs2,
#				 'figure.dpi':150,
#				 'figure.constrained_layout.use':True,
#				 'axes.edgecolor':'black',
#				 'axes.labelsize':fs1,
#				 'axes.labelcolor':'black',
#				 'axes.titlesize':fs2,
#				 'axes.titlepad':2,
#				 'axes.linewidth':0.5,
#				 'lines.linewidth':0.75,
#				 'lines.markersize': 2,
#				 'text.color':'black',
#				 'xtick.color':'black',
#				 'xtick.labelsize':fs1,
#				 'xtick.major.width':0.5,
#				 'xtick.major.size':3,
#				 'xtick.direction':'in',
#				 'ytick.color':'black',
#				 'ytick.labelsize':fs1,
#				 'ytick.major.width':0.5,
#				 'ytick.major.size':3,
#				 'ytick.direction':'in',
#				 'legend.fontsize':fs1,
#				 'savefig.dpi':900,
#				 'savefig.transparent':True,
#				 'savefig.format':'png',
#				 'savefig.pad_inches':0.1,
#				 'savefig.bbox':'tight'}
#	 elif type=='spyder_fs7':
#		 fs1=7
#		 fs2=7
#		 params={'font.sans-serif':'Arial',
#				 'font.size':fs1,
#				 'figure.titlesize':fs2,
#				 'figure.dpi':150,
#				 'figure.constrained_layout.use':True,
#				 'axes.edgecolor':'black',
#				 'axes.labelsize':fs1,
#				 'axes.labelcolor':'black',
#				 'axes.titlesize':fs2,
#				 'axes.titlepad':2,
#				 'axes.linewidth':0.5,
#				 'lines.linewidth':0.75,
#				 'lines.markersize': 2,
#				 'text.color':'black',
#				 'xtick.color':'black',
#				 'xtick.labelsize':fs1,
#				 'xtick.major.width':0.5,
#				 'xtick.major.size':3,
#				 'xtick.direction':'in',
#				 'ytick.color':'black',
#				 'ytick.labelsize':fs1,
#				 'ytick.major.width':0.5,
#				 'ytick.major.size':3,
#				 'ytick.direction':'in',
#				 'legend.fontsize':fs1,
#				 'savefig.dpi':900,
#				 'savefig.transparent':False,
#				 'savefig.format':'png',
#				 'savefig.pad_inches':0.1,
#				 'savefig.bbox':'tight'}
#	 elif type=='Presentation':
#		 fs1=7
#		 fs2=7
#		 params={'font.sans-serif':'Arial',
#				 'font.size':fs1,
#				 'figure.titlesize':fs2,
#				 'figure.dpi':150,
#				 'figure.constrained_layout.use':True,
#				 'axes.edgecolor':'black',
#				 'axes.labelsize':fs1,
#				 'axes.labelcolor':'black',
#				 'axes.titlesize':fs2,
#				 'axes.titlepad':2,
#				 'axes.linewidth':0.5,
#				 'lines.linewidth':1.25,
#				 'lines.markersize': 4,
#				 'text.color':'black',
#				 'xtick.color':'black',
#				 'xtick.labelsize':fs1,
#				 'xtick.major.width':0.5,
#				 'xtick.major.size':3,
#				 'xtick.direction':'in',
#				 'ytick.color':'black',
#				 'ytick.labelsize':fs1,
#				 'ytick.major.width':0.5,
#				 'ytick.major.size':3,
#				 'ytick.direction':'in',
#				 'legend.fontsize':fs1,
#				 'savefig.dpi':900,
#				 'savefig.transparent':True,
#				 'savefig.format':'png',
#				 'savefig.pad_inches':0.1,
#				 'savefig.bbox':'tight'}
#	 elif type=='Presentation Dark':

#		 fs1=7; fs2=7; cla=[0.92,0.92,0.92];
#		 params={'font.sans-serif':'Arial',
#				 'font.size':fs1,
#				 'figure.titlesize':fs2,
#				 'figure.dpi':150,
#				 'figure.constrained_layout.use':True,
#				 'axes.edgecolor':cla,
#				 'axes.labelsize':fs1,
#				 'axes.labelcolor':cla,
#				 'axes.titlesize':fs2,
#				 'axes.titlepad':2,
#				 'axes.linewidth':0.5,
#				 'lines.linewidth':5,
#				 'lines.markersize': 4,
#				 'text.color':cla,
#				 'xtick.color':cla,
#				 'xtick.labelsize':fs1,
#				 'xtick.major.width':0.5,
#				 'xtick.major.size':3,
#				 'xtick.direction':'in',
#				 'ytick.color':cla,
#				 'ytick.labelsize':fs1,
#				 'ytick.major.width':0.5,
#				 'ytick.major.size':3,
#				 'ytick.direction':'in',
#				 'legend.fontsize':fs1,
#				 'savefig.dpi':900,
#				 'savefig.transparent':True,
#				 'savefig.format':'png',
#				 'savefig.pad_inches':0.1,
#				 'savefig.bbox':'tight'}
#	 return params

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

		gp['printfigs']='On'
		gp['fig w mult']=1.0
		gp['fig h mult']=1.0
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
		gp['lw3']=1.5;
		gp['ms']=4
		gp['Alpha1']=0.225;
		gp['Alpha2']=0.45;
		gp['tickl']=1.5;

		params={'figure.dpi':125,
				'figure.titlesize':gp['fs2'],
				'figure.constrained_layout.use':False,
				'font.sans-serif':'Arial',
				'font.size':gp['fs1'],
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
				'savefig.bbox':None}

	elif type=='Web':

		#gp['printfigs']='Off'
		gp['fs_s']=5
		gp['fs_m']=6
		gp['fs_l']=7
		gp['fs_xl']=8
		gp['cla']=[0,0,0]
		gp['clt']=[0,0,0]
		gp['cl1']=[0.27,0.49,0.79];
		gp['cl2']=[0.6,0.9,0];
		gp['cl3']=[0.7,0.4,0.95];
		gp['cl4']=[0.85,0,0]
		gp['lw1']=0.5;
		gp['lw2']=1.0;
		gp['lw3']=1.5;
		gp['ms']=4
		gp['Alpha1']=0.225;
		gp['Alpha2']=0.45;
		gp['tickl']=1.5;

		params={'figure.dpi':125,
				'figure.constrained_layout.use':True,
				'figure.titlesize':gp['fs_m'],
				'font.size':gp['fs_m'],
				'xtick.labelsize':gp['fs_m'],
				'ytick.labelsize':gp['fs_m'],
				'axes.labelsize':gp['fs_m'],
				'axes.titlesize':gp['fs_l'],
				'legend.fontsize':gp['fs_m'],
				'font.sans-serif':'Arial',
				'axes.edgecolor':gp['cla'],
				'axes.labelcolor':gp['cla'],
				'axes.titlepad':2,
				'axes.linewidth':gp['lw1'],
				'lines.linewidth':gp['lw2'],
				'xtick.major.width':gp['lw1'],
				'lines.markersize':gp['ms'],
				'text.color':gp['cla'],
				'xtick.color':gp['cla'],
				'xtick.major.size':3,
				'xtick.direction':'in',
				'ytick.color':gp['cla'],
				'ytick.major.width':0.5,
				'ytick.major.size':3,
				'ytick.direction':'in',
				'savefig.dpi':900,
				'savefig.transparent':False,
				'savefig.format':'png',
				'savefig.pad_inches':0.1,
				'savefig.bbox':None}
	plt.rcParams.update(params)
	return gp

#%% Z-score a variable
# y,mu,sig=gu.zscore(x)
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

#%% Z-score a variable
# y,mu,sig=gu.zscore(x,ikp)
def zscore_ikp(x,ikp):
	if x.ndim==1:
		y=np.nan*np.ones(x.size)
		mu=np.nanmean(x[ikp])
		sig=np.nanstd(x[ikp])
		y[ikp]=(x[ikp]-mu)/sig
	else:
		y=np.nan*np.ones(x.shape)
		mu=np.zeros(x.shape[1])
		sig=np.zeros(x.shape[1])
		for j in range(x.shape[1]):
			mu[j]=np.nanmean(x[:,j])
			sig[j]=np.nanstd(x[:,j])
			y[:,j]=(x[:,j]-mu[j])/sig[j]
	return y,mu,sig

#%% Read NetCDF4
# nc=gu.ReadNC(fin)
def ReadNC(fin):
	ds=nc.Dataset(fin)
	d={}
	for k in ds.variables.keys():
		d[k]=np.array(ds.variables[k][:])
	return d

#%% Sum over interval
def BlockSum(x,ivl):
	if x.ndim==1:
		y=np.zeros(int(x.size/ivl))
		cnt=0
		for i in range(0,x.size,ivl):
			y[cnt]=np.sum(x[i:i+ivl])
			cnt=cnt+1
	elif x.ndim==2:
		y=np.zeros((int(x.shape[0]/ivl),x.shape[1]))
		cnt=0
		for i in range(0,x.shape[0],ivl):
			y[cnt,:]=np.sum(x[i:i+ivl,:],axis=0)
			cnt=cnt+1
	else:
		y=np.zeros((int(x.shape[0]/ivl),x.shape[1],x.shape[2]))
		cnt=0
		for i in range(0,x.shape[0],ivl):
			y[cnt,:,:]=np.sum(x[i:i+ivl,:,:],axis=0)
			cnt=cnt+1
	return y

#%% Sum over interval
def BlockMean(x,ivl):
	if x.ndim==1:
		c=np.arange(0,x.size,ivl,dtype=int)
		y=np.zeros(c.size)
		for i in range(c.size):
			y[i]=np.mean(x[c[i]:c[i]+ivl])
	elif x.ndim==2:
		y=np.zeros((int(x.shape[0]/ivl),x.shape[1]))
		cnt=0
		for i in range(0,x.shape[0],ivl):
			y[cnt,:]=np.mean(x[i:i+ivl,:],axis=0)
			cnt=cnt+1
	else:
		y=np.zeros((int(x.shape[0]/ivl),x.shape[1],x.shape[2]))
		cnt=0
		for i in range(0,x.shape[0],ivl):
			y[cnt,:,:]=np.mean(x[i:i+ivl,:,:],axis=0)
			cnt=cnt+1
	return y

#%% Sum over interval
def BlockMissing(x,ivl):
	c=np.arange(0,x.size,ivl,dtype=int)
	y=np.zeros(c.size)
	for i in range(c.size):
		x0=x[c[i]:c[i]+ivl]
		ind=np.where(np.isnan(x0)==True)[0]
		y[i]=ind.size
	return y

#%% Inverval Mean
def IvlMean(x,ivl):
	if x.ndim==1:
		y=np.zeros(ivl)
		for i in range(ivl):
			y[i]=np.nanmean(x[i::ivl])
	elif x.ndim==2:
		y=np.zeros((ivl,x.shape[1]))
		for i in range(ivl):
			y[i,:]=np.nanmean(x[i::ivl,:],axis=0)
	else:
		y=[]
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
	x=np.array(x)
	if x.size==1:
		y=np.maximum(np.minimum(xmax,x),xmin)
	else:
		y=np.maximum(np.minimum(xmax,x),xmin)

	return y

#%% Density profile
def ksdensity(y):
	kde=stats.gaussian_kde(y)
	x=np.linspace(y.min(),y.max(),100)
	p=kde(x)
	return p

#%% Polynomial Surface Fit
def PolynomialSurfaceFit(x,y,z,xI,yI,maskI,beta):
	
	#x=x0/1000
	#y=y0/1000
	#z=z0
	
	# Parameters
	m,n=xI.shape
	L=z.size
	zI=np.zeros(shape=(m,n))
	eI=np.zeros(shape=(m,n))
	nI=np.zeros(shape=(m,n))
		
	nth=40 # Maximum number of neighbouring measurements used for interpolation
	#tol=[1e-12]; % Tolerance for the rejection of small singular values
	tol=1e-15
	sfar=15000 # Maximum separation
	
	def fhyp(b,x):
		y=((b[0]*b[1]*x)/(b[0]*x+b[2]))+b[2]
		return y
	flg=0
	if flg==1:
		beta=np.array([0.00001,0.25,0.05])
		d0=np.arange(0,100)
		s=fhyp(beta,d0*1000)
		plt.close('all'); plt.plot(d0,s,'b-')

	for i in range(m):
		for j in range(n):
			if maskI[i,j]==0:
				continue
			
			# Calculate distance
			d0=np.sqrt( (xI[i,j]*np.ones(L)-x)**2 + (yI[i,j]*np.ones(L)-y)**2 )
			
			# Put them in order of increasing distance
			ord=np.argsort(d0)
			ord=ord[1:np.minimum(nth,ord.size)]
			d0=d0[ord]
			x0=x[ord]
			y0=y[ord]
			z0=z[ord]
			#plt.plot(d0)
			
			# Exclude points too far away
			iKeep=np.where(d0<sfar)[0]
			if iKeep.size<12:
				continue
			d0=d0[iKeep]
			x0=x0[iKeep]
			y0=y0[iKeep]
			z0=z0[iKeep]
		
			# Calculate weights based on errorgram (standard error versus separation)
			s=fhyp(beta,d0*1000)
			#s=0.0001*np.ones(z0.size)
			#print(np.mean(d0))
			
			# Interpolate
			b,e,istat=polyfit_lls(z0,s,x0-xI[i,j],y0-yI[i,j],tol)
 
			zI[i,j]=b[0] 
			#eI[i,j]=e[0]
			nI[i,j]=iKeep.size

	return zI,eI,nI

#%% Linear least squares polynomial fit based on singular value decomposition
def polyfit_lls(z,s,x,y,tol):
	# istat is an output status code:
	#   0 normal completion
	#   1 M > N
	#   2 a standard error S was found <=0
	#   3 error using SVD
	istat=0
	n=z.size
	
	# Create design matrix BAS and the number of basis functions M based on L
	# which is determined by the number of available neighbouring measurements
	if n>=28*2:
		l=6
	elif n>=21*2:
		l=5
	elif n>=15*2:
		l=4
	elif n>=10*2:
		l=3
	elif n>=6*2:
		l=2
	else:
		# n <= m
		istat=1
		return
		
	m,bas=BasisN(x,y,l)
	
	# Weighted Z ZETA and weighted average of Z ZETAMU
	s2=(1/s)**2
	zeta=(1/s)*z
	zetamu=np.sum(s2*z)/s2
	
	# Weighted design matrix DES	
	a=np.tile(1/s,(m,1)).T
	des=a*bas
	
	if np.sum(np.isnan(des))>0:
		print('Stopping')
		return

	# Estimate parameers by singular value decomposition
	u,w,v=np.linalg.svd(des)
	
	# Scan for singluar values < TOL
	#w=np.diag(w)
	iToss=np.where(w<tol*np.max(w))
	w[iToss]=0
	
	# Solve U*B=ZETA for the parameter vector B, where U has been decomposed 
	# into U,W, and V.	
	work=np.nan*np.ones(w.size)
	b=np.nan*np.ones(w.size)   
	for j in range(m):
		work[j]=0
		if w[j]!=0:
			for i in range(n):
				work[j]=work[j]+u[i,j]*zeta[i]					
			work[j]=work[j]/w[j]
	vT=v.T
	for j in range(m):
		b[j]=0
		for k in range(m):
			b[j]=b[j]+vT[j,k]*work[k]
	
	# Goodness of fit of the polynomial at the points in the scatter
	dz=(z-zetamu)/s
	vobs=np.sum(dz**2)
	bt=np.tile(b,(n,1))
	zhat=np.sum(bt*bas,axis=1)
	dz=(z-zhat)/s
	vfit=np.sum(dz**2)
	vfit=100*(1-vfit/vobs)
	   
	# Calculate covariance matrix C and standard errors E of the parameters B
	#wf=np.zeros((w.size,w.size))
	#wf=np.fill_diagonal(wf,w)
	#c,e=svdcovm(wf,v,m)
	e=1
	
	return b,e,istat

#%% Covariance matrix and parameter standard errors
def svdcovm(w,v,m):
	c=np.nan*np.ones((m,m))
	e=np.nan*np.ones((m,m))
	ind=np.where(w!=0)
	w[ind]=1/(w**2)
	for i in range(m):
		for j in range(i):
			s=0
			for k in range(m):
				s=s+v[i,k]*v[j,k]*w[k]
			c[i,j]=s
			c[j,i]=s	
	for i in range(m):
		if c[i,i]<0:
			e[i]=999
		else:
			e[i]=np.sqrt(c[i,i])
	return c,e

#%% Basis functions for bivariate polynomials
def BasisN(x,y,n):	
	# Column dimension of design matrix
	l=int((n+1)*(n+2)/2)
	# Design matrix
	g0=np.ones(x.size)	
	if n==0:
		g=g0
	elif n==1:
		g=[g0,x,y]
	elif n==2:
		g=[g0,x,y,x**2,x*y,y**2]
	elif n==3:
		g=[g0,x,y,x**2,x*y,y**2,x**3,x*2*y,x*y**2,y**3]
	elif n==4:
		g=[g0,x,y,x**2,x*y,y**2,x**3,x**2*y,x*y**2,y**3,x**4,x**3*y,x**2*y**2,x*y**3,y**4]
	elif n==5:
		g=[g0,x,y,x**2,x*y,y**2,x**3,x**2*y,x*y**2,y**3,x**4,x**3*y,x**2*y**2,x*y**3,y**4,x**5,x**4*y,x**3*y**2,x**2*y**3,x*y**4,y**5]
	elif n==6:
		g=[g0,x,y,x**2,x*y,y**2,x**3,x**2*y,x*y**2,y**3,x**4,x**3*y,x**2*y**2,x*y**3,y**4,x**5,x**4*y,x**3*y**2,x**2*y**3,x*y**4,y**5,x**6,x**5*y,x**4*y**2,x**3*y**3,x**2*y**4,x*y**5,y**6]
	g=np.array(g).T
	return l,g

#%% Bunch contents of dictionary
class Bunch(object):
	def __init__(self, adict):
		self.__dict__.update(adict)

#%%
# z=gu.GeneratePerlin(shape,p_th)
def GeneratePerlin(shape,p_th):	
	x=np.linspace(0,10,shape[1])
	y=np.linspace(0,10,shape[0])
	xv,yv=np.meshgrid(x,y)
	seed=np.random.randint(1,500,1)[0]
	noise=perlin(xv,yv,seed=seed)	
	p=np.percentile(noise,p_th)
	z=np.zeros(noise.shape,dtype='int8')
	ind=np.where(noise<p); z[ind]=1
	#plt.close('all'); plt.matshow(z)
	return z
def perlin(x,y,seed=0):
	np.random.seed(seed)
	p = np.arange(256, dtype=int)
	np.random.shuffle(p)
	p = np.stack([p, p]).flatten()
	xi = x.astype(int)
	yi = y.astype(int)
	xf = x - xi
	yf = y - yi
	u = fade(xf)
	v = fade(yf)
	n00 = gradient(p[p[xi] + yi], xf, yf)
	n01 = gradient(p[p[xi] + yi + 1], xf, yf - 1)
	n11 = gradient(p[p[xi + 1] + yi + 1], xf - 1, yf - 1)
	n10 = gradient(p[p[xi + 1] + yi], xf - 1, yf)
	x1 = lerp(n00, n10, u)
	x2 = lerp(n01, n11, u)
	return lerp(x1, x2, v)
def fade(t):
	return 6 * t**5 - 15 * t**4 + 10 * t**3
def lerp(a, b, x):
	return a + x * (b - a)
def gradient(h, x, y):
	vectors = np.array([[0, 1], [0, -1], [1, 0], [-1, 0]])
	g = vectors[h % 4]
	return g[:,:,0] * x + g[:,:,1] * y

#%%
def ScatterDensity(x,y,limX,limY,nBx,nBy):
	x=x.flatten()
	y=y.flatten()
	binX=np.linspace(limX[0],limX[1],nBx)
	binY=np.linspace(limY[0],limY[1],nBy)
	dX=binX[1]-binX[0]
	dY=binY[1]-binY[0]
	Z=np.zeros((binY.size,binX.size))
	for i in range(binY.size):
		for j in range(binX.size):
			ind=np.where( (np.abs(x-binX[i])<dX/2) & (np.abs(y-binY[j])<dY/2) )[0]
			Z[i,j]=ind.size
	return binX,binY,Z

#%% Get mean values of y for unique values of x
def mean_ucat(xcat,y):
	u=np.unique(xcat)
	yU=np.zeros(u.size)
	for i in range(u.size):
		ind=np.where(xcat==u[i])[0]
		yU[i]=np.mean(y[ind])
	return u,yU

#%%
def Month2DOY(month):
	dim=np.array([31,28,31,30,31,30,31,31,30,31,30,31])
	cdim=np.cumsum(dim)
	doy=cdim[month-1]
	return doy

#%% Used to unpack bits from floating point Geotiff from MODIS USGS
def unpackbits(x, num_bits):
	if np.issubdtype(x.dtype, np.floating):
		raise ValueError("numpy data type needs to be int-like")
	xshape = list(x.shape)
	x = x.reshape([-1, 1])
	mask = 2**np.arange(num_bits, dtype=x.dtype).reshape([1, num_bits])
	return (x & mask).astype(bool).astype(int).reshape(xshape + [num_bits])

#%%