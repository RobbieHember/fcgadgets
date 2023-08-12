
import numpy as np
from fcgadgets.cbrunner import cbrun_util as cbu

#%% Schedule non-obligation stand establishment

def PredictNOSE_OnTheFly(meta,pNam,iScn,iBat,vi,iT):

    if meta[pNam]['Year'][iT]<meta[pNam]['Project']['Year End']-3:
        
        # Initialize project type
        if 'RegTypeNO' not in meta[pNam]['Project'].keys():
            meta[pNam]['Project']['RegTypeNO']=np.zeros(meta[pNam]['Project']['N Stand'],dtype='int16')
    
        flg=np.zeros((meta[pNam]['Project']['Batch Size'][iBat],meta['Core']['Max Events Per Year']))
        ind=np.where(vi['EC']['ID Event Type'][iT-2,:,:]==meta['LUT']['Event']['Wildfire'])
        flg[ind]=1
        flg=np.max(flg,axis=1)
        #print(flg.shape)
    
        Po=0.03 #meta[pNam]['Scenario'][iScn]['NOSE Prob']
        rn=np.random.random(flg.size)   
    
        indAffected=np.where( (rn<Po) & (flg==1) )[0]
        #print(indAffected.size)
        if indAffected.size>0:
            for iA in indAffected:
                
                # Add wildfire at 100% mortality to both scenarios
                iAvailable=np.where(vi['EC']['ID Event Type'][iT,iA,:]==0)[0]
                if iAvailable.size>0:
                    iE=iAvailable[0]
                    vi['EC']['ID Event Type'][iT,iA,iE]=meta['LUT']['Event']['Wildfire']
                    vi['EC']['Mortality Factor'][iT,iA,iE]=1.0
                    vi['EC']['ID Growth Curve'][iT,iA,iE]=1
                
                # Add planting to actual scenario, natural regen to baseline to next year
                iT2=iT+1
                iAvailable=np.where(vi['EC']['ID Event Type'][iT2,iA,:]==0)[0]
                if iAvailable.size>0:
                    iE=iAvailable[0]
                    meta[pNam]['Project']['RegTypeNO'][ meta[pNam]['Project']['indBat'][iA] ]=meta['LUT']['Derived']['RegenTypeNO']['Straight Fire']
                    meta[pNam]['Project']['Strata']['Project Type']['ID'][ meta[pNam]['Project']['indBat'][iA] ]=meta['LUT']['Derived']['RegenTypeNO']['Straight Fire']
                    if np.isin(iScn,meta[pNam]['Project']['Baseline Indices'])==True:                  
                        vi['EC']['ID Event Type'][iT2,iA,iE]=meta['LUT']['Event']['Regen at 25% Growth']
                        vi['EC']['Mortality Factor'][iT2,iA,iE]=0
                        #vi['EC']['Growth Factor'][iT2,iA,iE]=75
                        vi['EC']['ID Growth Curve'][iT2,iA,iE]=1
                    else:                       
                        vi['EC']['ID Event Type'][iT2,iA,iE]=meta['LUT']['Event']['Planting']
                        vi['EC']['Mortality Factor'][iT2,iA,iE]=1
                        vi['EC']['ID Growth Curve'][iT2,iA,iE]=2
                else:
                    print('No space left in event chronology for on-the-fly breakup!')

    return vi    