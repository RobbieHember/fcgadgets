'''
UNCERTAINTY FRAMEWORK
'''

import numpy as np

def normal(x,mu,sigma):
	return ( 2.*np.pi*sigma**2. )**-.5 * np.exp( -.5 * (x-mu)**2. / sigma**2. )

mu,sigma,n=-20,20,1000
x=np.random.normal(mu,sigma,n) #generate random list of points from normal distribution
y=normal(x,mu,sigma)

ord=np.argsort(x)
x=x[ord]
y=y[ord]

#%% Plot distribution

cl1=[0.24,0.49,0.76]
cl2=[0.8,0.45,0.2]

plt.close('all'); fig,ax=plt.subplots(1,figsize=gu.cm2inch(15,8))
ax.plot(x,y,color=cl1)
ax.plot([mu,mu],[0,0.03],'k--',color=cl1)
ax.plot([0,0],[0,0.03],'k--',color=cl2)
ax.set(position=[0.08,0.15,0.84,0.8],ylabel='Frequency',xlabel='Response of heterotrophic respiration\nto N application (%)',ylim=[0,0.025])
ax.yaxis.set_ticks_position('both'); ax.xaxis.set_ticks_position('both'); ax.tick_params(length=meta['Graphics']['gp']['tickl'])
ax.text(-20,0.01,'Review of\n literature',fontsize=6,ha='center',backgroundcolor='w',color=cl1)
ax.text(0,0.019,'Exclusion of\nthe process\n(conservativeness\nprinciple)',fontsize=6,ha='center',backgroundcolor='w',color=cl2)
#plt.axis_tight()
gu.PrintFig(r'C:\Users\rhember\OneDrive - Government of BC\Figures\Mitigation\UncertaintyDistributionExample','png',900)

#%%

