'''
Markov Chain Monte Carlo (MCMC) Utilities
'''

#%% Import modules
import numpy as np
import scipy.stats as sc
from matplotlib import pyplot as plt
from tqdm import tqdm
import seaborn as sns

#%% Calculate log likelihood
# Assumes the last parameter is the error variance
def log_likelihood(theta,y,yhat):    
    ll=np.sum(sc.norm.logpdf(y,loc=yhat,scale=theta[-1]))
    return ll

#%% Calculate log prior
def log_prior(theta,thetaScale):
    l_prior=0
    for i in range(theta.size):
        l_prior=l_prior+sc.norm.logpdf(theta[i],loc=0,scale=thetaScale[i])
    return l_prior

#%% Proposed new parameters
def proposal(thetaPrv,SearchWidth,thetaMin,thetaMax):
    thetaProposed=sc.norm.rvs(loc=thetaPrv,scale=SearchWidth,size=(thetaPrv.shape))
    thetaProposed=np.maximum(thetaMin,thetaProposed)
    thetaProposed=np.minimum(thetaMax,thetaProposed)
    return thetaProposed

#%% Metropolis Hastings version of MCMC (with yobs input to function)
def MetropolisHastings(thetaInit,thetaScale,thetaMin,thetaMax,SearchWidth,y,model,N_Iter):
    chain=np.zeros((N_Iter,thetaInit.size))
    chain[0]=thetaInit
    ll=np.zeros(N_Iter)
    ar=np.zeros(N_Iter)
    for i in tqdm(range(N_Iter-1)):        
        thetaProposed=proposal(chain[i],SearchWidth,thetaMin,thetaMax)        
        yhatProposed,Dummy=model(thetaProposed)        
        log_post_proposal=log_prior(thetaProposed,thetaScale)+log_likelihood(thetaProposed,y,yhatProposed)        
        if i==0:
            yhatPrv,Dummy=model(chain[i])
            log_post_prv=log_prior(chain[i],thetaScale)+log_likelihood(chain[i],y,yhatPrv)        
        post_prob=np.exp(log_post_proposal-log_post_prv) # symmetric proposal
        if np.random.rand() < post_prob:
            chain[i+1]=thetaProposed
            ll[i+1]=log_post_proposal
            log_post_prv=log_post_proposal
            yhatPrv=yhatProposed
            ar[i]=1
        else:
            chain[i+1]=chain[i]
            ll[i+1]=log_post_prv        
    return chain,ll,ar

#%% Metropolis Hastings version of MCMC (with model outputting yobs)
def MetropolisHastings2(thetaInit,thetaScale,thetaMin,thetaMax,SearchWidth,y,model,N_Iter):
    chain=np.zeros((N_Iter,thetaInit.size))
    chain[0]=thetaInit
    ll=np.zeros(N_Iter)
    ar=np.zeros(N_Iter)
    for i in tqdm(range(N_Iter-1)):        
        thetaProposed=proposal(chain[i],SearchWidth,thetaMin,thetaMax)        
        yhatProposed,y=model(thetaProposed)        
        log_post_proposal=log_prior(thetaProposed,thetaScale)+log_likelihood(thetaProposed,y,yhatProposed)        
        if i==0:
            yhatPrv,y=model(chain[i])
            log_post_prv=log_prior(chain[i],thetaScale)+log_likelihood(chain[i],y,yhatPrv)        
        post_prob=np.exp(log_post_proposal-log_post_prv) # symmetric proposal
        if np.random.rand()<post_prob:
            chain[i+1]=thetaProposed
            ll[i+1]=log_post_proposal
            log_post_prv=log_post_proposal
            yhatPrv=yhatProposed
            ar[i]=1
        else:
            chain[i+1]=chain[i]
            ll[i+1]=log_post_prv
    return chain,ll,ar

#%% Metropolis Hastings version of MCMC (with model outputting yobs,potentially not starting from scratch)
def MetropolisHastings2b(thetaInit,thetaScale,thetaMin,thetaMax,SearchWidth,y,model,N_Iter,chain,ll,ar,E,iStart):    
    for i in tqdm(range(iStart,N_Iter-1)):        
        thetaProposed=proposal(chain[i],SearchWidth,thetaMin,thetaMax)        
        yhatProposed,y,E[i,:]=model(thetaProposed)        
        log_post_proposal=log_prior(thetaProposed,thetaScale)+log_likelihood(thetaProposed,y,yhatProposed)        
        if i==iStart:
            yhatPrv,y,E[i]=model(chain[i])
            log_post_prv=log_prior(chain[i],thetaScale)+log_likelihood(chain[i],y,yhatPrv)        
        post_prob=np.exp(log_post_proposal - log_post_prv)  # symmetric proposal
        if np.random.rand() < post_prob:
            chain[i+1]=thetaProposed
            ll[i+1]=log_post_proposal
            log_post_prv=log_post_proposal
            yhatPrv=yhatProposed
            ar[i]=1
        else:
            chain[i+1]=chain[i]
            ll[i+1]=log_post_prv
        print(ll[i+1])
    return chain,ll,ar,E

#%% Test

def Test():
    
    global x,y,yhat
    
    # Training data
    thetaTrue=[0,5,10]
    x=np.linspace(-20,20)    
    y=thetaTrue[0]+thetaTrue[1]*x+sc.norm.rvs(loc=0,scale=thetaTrue[-1],size=x.shape)
    #plt.close('all'); plt.plot(x,y,'.'); plt.xlabel('x'); plt.ylabel('y'); plt.title('data')
    
    # our linear model a+b*x
    def model(theta):        
        yhat=theta[0]+theta[1]*x
        Dummy=[]
        return yhat,Dummy

    #thetaInit=thetaTrue
    thetaInit=np.array([4,10,5])
    thetaScale=np.array([5,10,30])
    thetaMin=thetaInit-5*thetaInit
    thetaMax=thetaInit+5*thetaInit
    SearchWidth=np.array([0.25,0.5,1.0])
    
    N_Iter=10000
    
    chain,ll,ar=MetropolisHastings(thetaInit,thetaScale,thetaMin,thetaMax,SearchWidth,y,model,N_Iter)

    # Analyze results
    plt.close('all')
    fig,ax=plt.subplots(1)
    plt.plot(ll,'-o')

    burn_in=int(0.5*N_Iter)
    fig,ax=plt.subplots(1,3)
    for i in range(len(ax)):
        sns.distplot(chain[burn_in:,i],ax=ax[i])
        ax[i].set_xlabel(f"param[{i}]")
    ax[0].set_ylabel('posterior density')
    fig.tight_layout()

    return