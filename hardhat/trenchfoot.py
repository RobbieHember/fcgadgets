
import numpy as np
from fcgadgets.cbrunner import cbrun_util as cbu

#%% Schedule non-obligation stand establishment

def PredictNOSE_OnTheFly(meta,pNam,iScn,iBat,vi,iT):

    if meta[pNam]['Year'][iT]<meta[pNam]['Project']['Year End']-3:
        
        # Initialize project type
        if 'RegType' not in meta[pNam]['Project'].keys():
            meta[pNam]['Project']['RegType']=np.zeros(meta[pNam]['Project']['N Stand'],dtype='int16')
    
        # Identify wildfire occurrence
        flg_wf=np.zeros((meta[pNam]['Project']['Batch Size'][iBat],meta['Core']['Max Events Per Year']))
        ind=np.where(vi['EC']['ID Event Type'][iT-2,:,:]==meta['LUT']['Event']['Wildfire'])
        flg_wf[ind]=1
        flg_wf=np.max(flg_wf,axis=1)
        #print(flg_wf.shape)
    
        Po=meta[pNam]['Scenario'][iScn]['NOSE Prob']
        rn=np.random.random(flg_wf.size)   
    
        indAffected=np.where( (rn<Po) & (flg_wf==1) )[0]
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
                    meta[pNam]['Project']['RegType'][ meta[pNam]['Project']['indBat'][iA] ]=meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Wildfire']
                    meta[pNam]['Project']['Strata']['Project Type']['ID'][ meta[pNam]['Project']['indBat'][iA] ]=meta['LUT']['Derived']['RegenType']['Straight-to-planting Post Wildfire']
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