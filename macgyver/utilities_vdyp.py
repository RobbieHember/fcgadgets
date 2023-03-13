
import numpy as np
import pandas as pd
import pylab as pl
import matplotlib.pyplot as plt
import time
import fcgadgets.macgyver.utilities_general as gu
import fcexplore.psp.Processing.psp_utilities as utl_gp

#%% Set figure properties

gp=gu.SetGraphics('Manuscript')

#%% Impot ground plots

meta={}
meta['Paths']={}
meta['Paths']['DB']=r'C:\Users\rhember\Documents\Data\GroundPlots\PSP-NADB2'
meta=utl_gp.ImportParameters(meta)
d=gu.ipickle(meta['Paths']['DB'] + '\\Processed\\L2\\L2_BC.pkl')
sobs=d['sobs'].copy(); del d

#%% Prepare VDYP inputs from ground plots

def PrepareInputsFromGroundPlots():

    # Notes:
    # 1) Problems running model? Confirm that schema file is in the input folder

    pthout=r'D:\Data\FCI_Projects\BCFCS_EvalAtPlots\Inputs\VDYP'

    m=sobs['Year t0'].size
    id=np.arange(1,m+1,1)

    # Polygon ID
    nams=['POLYGON_RCRD_ID','MAINTAINER','MAP_ID','MAP_QUAD','MAP_SUB_QUAD','POLYGON_ID']
    df1=pd.DataFrame(np.zeros(shape=(m,len(nams))),columns=nams)
    df1.iloc[:,0]=id                    # Polygon ID or df0['POLYGON_ID']
    df1.iloc[:,1]=1                     # Maintainer, district string e.e. DPG
    df1.iloc[:,2]=''                    # Map id
    df1.iloc[:,3]=0                     # Not used
    df1.iloc[:,4]=0                     # Not used
    df1.iloc[:,5]=id                    # Polygon ID
    df1.to_csv(pthout + '\\POLYGON_ID.csv',sep=',',index=False)

    # Polygon
    nams=['POLYGON_RCRD_ID','INVENTORY_STANDARD','REFERENCE_YEAR','FIZ','BEC_ZONE','PCT_STOCKABLE','NON_PROD_DESC','YIELD_FACTOR']
    df1=pd.DataFrame(np.zeros(shape=(m,len(nams))),columns=nams)
    df1.iloc[:,0]=id                         # Polygon ID
    df1.iloc[:,1]='F'                        # Inventory standar
    df1.iloc[:,2]=sobs['Year t0'] # Measurement (reference) year ???
    df1.iloc[:,3]=sobs['FIZ'] # FIZ
    df1.iloc[:,4]=utl_gp.lut_id2cd(meta,'Ecozone BC L1',sobs['Ecozone BC L1']) # BEC zone
    df1.iloc[:,5]=''                         # Stockable normally not supplied
    df1.iloc[:,6]=''                         # Non prod desc
    df1.iloc[:,7]=''                         # Yield factor normally not supplied
    df1.to_csv(pthout + '\\POLYGON.csv',sep=',',index=False)

    # Layer
    nams=['LAYER_RCRD_ID','POLYGON_RCRD_ID','LAYER_ID','VDYP7_LAYER_ID','RANK_CODE','NON_FOREST_DESC','MEAS_UTIL_LEVEL','CC','BA','TPH','EST_SI_SPCS','EST_SI']
    df1=pd.DataFrame(np.zeros(shape=(m,len(nams))),columns=nams)
    df1.iloc[:,0]=id                                # Layer RCRD ID
    df1.iloc[:,1]=id                                # Polygon ID
    df1.iloc[:,2]=1                                 # Layer ID
    df1.iloc[:,3]='P'                               #
    df1.iloc[:,4]=1                                 # Rank code
    df1.iloc[:,5]=''                                #
    df1.iloc[:,6]=7.5                               # Not used, model will not run if you enter anything but 7.5, fill with 7.5
    df1.iloc[:,7]=''                                # Crown closure
    df1.iloc[:,8]=sobs['BA L t0']                   # Basal area
    df1.iloc[:,9]=sobs['N L t0']      # Live stems per hectare
    df1.iloc[:,10]='' #df0['EST_SITE_INDEX_SPECIES_CD'] # Estimated SI species
    df1.iloc[:,11]='' #df0['SITE_INDEX'] # Estimated SI, exclude in BAH???
    df1.to_csv(pthout + '\\LAYER.csv',sep=',',index=False)

    # Species
    nams=['SPECIES_RCRD_ID','LAYER_RCRD_ID','SPECIES_ID','SPECIES_CD','SPECIES_PCT','AGE','HEIGHT','SI','YTBH','BYAGE','SITE_CURVE']
    df1=pd.DataFrame(np.zeros(shape=(m*4,len(nams))),columns=nams)
    cnt=0
    for i in range(len(df0)):
        for j in range(4):
            tmp=np.isnan(df0['SPECIES_PCT_' + str(j+1)][i])
            if tmp==True:
                continue
            df1.iloc[cnt,0]=id[i]                              # Species rcrd id
            df1.iloc[cnt,1]=id[i]                              # layer rcrd id
            df1.iloc[cnt,2]=j+1                                # species id
            df1.iloc[cnt,3]=df0['SPECIES_CD_' + str(j+1)][i]   # species cd
            df1.iloc[cnt,4]=df0['SPECIES_PCT_' + str(j+1)][i]  # species pct
            df1.iloc[cnt,5]=df0['PROJ_AGE_1'][i]               # Age
            df1.iloc[cnt,6]=''                                 # Height %h_sps; % (include in BAH)
            df1.iloc[cnt,7]=df0['SITE_INDEX'][i]               # SI
            df1.iloc[cnt,8]=10                                 # ytbh
            df1.iloc[cnt,9]=12                                 # bhage
            df1.iloc[cnt,10]=''                                # site curve
            cnt=cnt+1
    df1.to_csv(pthout + '\\SPECIES.csv',sep=',',index=False)

    # Non-veg
    nams=['POLYGON_RCRD_ID','NON_VEG_ID','NON_VEG_COVER_TYPE','NON_VEG_COVER_PCT']
    df1=pd.DataFrame(np.zeros(shape=(m,len(nams))),columns=nams)
    df1.iloc[:,0]=id                                # polygon id
    df1.iloc[:,1]=id                                # non veg id
    df1.iloc[:,2]='BR'                              # non veg type
    df1.iloc[:,3]=0                                 # non veg percent
    df1.to_csv(pthout + '\\NON_VEG.csv',sep=',',index=False)

    # History
    nams=['HISTORY_RCRD_ID','LAYER_RCRD_ID','HISTORY_ID','SILV_BASE_CD','START_YEAR','END_YEAR','LAYER_PCT']
    df1=pd.DataFrame(np.zeros(shape=(m,len(nams))),columns=nams)
    df1.iloc[:,0]=id                                # History rcrd id
    df1.iloc[:,1]=id                                # Layer rcrd id
    df1.iloc[:,2]=id                                # History id
    df1.iloc[:,3]=''                                # SILV Base
    df1.iloc[:,4]=''                                # Start year ?
    df1.iloc[:,5]=''                                # End year ?
    df1.iloc[:,6]=''                                # Layer percent ?
    df1.to_csv(pthout + '\\HISTORY.csv',sep=',',index=False)

    # VRI Adjust
    nams=['VRIADJST_RCRD_ID','POLYGON_RCRD_ID','LAYER_ID','LH_075','BA_125','WSV_075','WSV_125','VCU_125','VD_125','VDW_125']
    df1=pd.DataFrame(np.zeros(shape=(m,len(nams))),columns=nams)
    df1.iloc[:,0]=''                 # vriadj rcrd id
    df1.iloc[:,1]=''                 # polygon rcrd id
    df1.iloc[:,2]=1                                 # layer id
    df1.iloc[:,3]=''                                #
    df1.iloc[:,4]=''                                #
    df1.iloc[:,5]=''                                #
    df1.iloc[:,6]=''                                #
    df1.iloc[:,7]=''                                #
    df1.iloc[:,8]=''                                #
    df1.iloc[:,9]=''                                #
    df1.to_csv(pthout + '\\VRIADJST.csv',sep=',',index=False)

    return



#%%

pthin=r'C:\Users\rhember\Documents\Data\GHG_DST\FCI_DST_RehabFire\GIS\ForTIPSYandCBM3.xlsx'
pthout=r'C:\Users\rhember\Documents\Data\GHG_DST\FCI_DST_RehabFire\VDYP\Inputs'

def PrepareInputsForVDYP(pthin,pthout):

    # Notes:
    # 1) Problems running model? Confirm that schema file is in the input folder

    # Import pre-fire stand info
    df0=pd.read_excel(pthin)
    m=df0.shape[0]
    id=np.arange(1,m+1,1)

    # Polygon ID
    nams=['POLYGON_RCRD_ID','MAINTAINER','MAP_ID','MAP_QUAD','MAP_SUB_QUAD','POLYGON_ID']
    df1=pd.DataFrame(np.zeros(shape=(len(df0),len(nams))),columns=nams)
    df1.iloc[:,0]=id                    # Polygon ID or df0['POLYGON_ID']
    df1.iloc[:,1]=1                     # Maintainer, district string e.e. DPG
    df1.iloc[:,2]=''                    # Map id
    df1.iloc[:,3]=0                     # Not used
    df1.iloc[:,4]=0                     # Not used
    df1.iloc[:,5]=id                    # Polygon ID
    df1.to_csv(pthout + '\\POLYGON_ID.csv',sep=',',index=False)

    # Polygon
    nams=['POLYGON_RCRD_ID','INVENTORY_STANDARD','REFERENCE_YEAR','FIZ','BEC_ZONE','PCT_STOCKABLE','NON_PROD_DESC','YIELD_FACTOR']
    df1=pd.DataFrame(np.zeros(shape=(len(df0),len(nams))),columns=nams)
    df1.iloc[:,0]=id                         # Polygon ID
    df1.iloc[:,1]='F'                        # Inventory standar
    df1.iloc[:,2]=df0['REFERENCE_YEAR']      # Measurement (reference) year ???
    df1.iloc[:,3]=df0['FIZ_CD']              # FIZ
    df1.iloc[:,4]=df0['BEC_ZONE_CODE']       # BEC zone
    df1.iloc[:,5]=''                         # Stockable normally not supplied
    df1.iloc[:,6]=''                         # Non prod desc
    df1.iloc[:,7]=''                         # Yield factor normally not supplied
    df1.to_csv(pthout + '\\POLYGON.csv',sep=',',index=False)

    # Layer
    nams=['LAYER_RCRD_ID','POLYGON_RCRD_ID','LAYER_ID','VDYP7_LAYER_ID','RANK_CODE','NON_FOREST_DESC','MEAS_UTIL_LEVEL','CC','BA','TPH','EST_SI_SPCS','EST_SI']
    df1=pd.DataFrame(np.zeros(shape=(len(df0),len(nams))),columns=nams)
    df1.iloc[:,0]=id                                # Layer RCRD ID
    df1.iloc[:,1]=id                                # Polygon ID
    df1.iloc[:,2]=1                                 # Layer ID
    df1.iloc[:,3]='P'                               #
    df1.iloc[:,4]=1                                 # Rank code
    df1.iloc[:,5]=''                                #
    df1.iloc[:,6]=7.5                               # Not used, model will not run if you enter anything but 7.5, fill with 7.5
    df1.iloc[:,7]=''                                # Crown closure
    df1.iloc[:,8]=''                                # Basal area
    df1.iloc[:,9]=df0['VRI_LIVE_STEMS_PER_HA']      # Live stems per hectare
    df1.iloc[:,10]=df0['EST_SITE_INDEX_SPECIES_CD'] # Estimated SI species
    df1.iloc[:,11]=df0['SITE_INDEX']                # Estimated SI, exclude in BAH???
    df1.to_csv(pthout + '\\LAYER.csv',sep=',',index=False)

    # Species
    nams=['SPECIES_RCRD_ID','LAYER_RCRD_ID','SPECIES_ID','SPECIES_CD','SPECIES_PCT','AGE','HEIGHT','SI','YTBH','BYAGE','SITE_CURVE']
    df1=pd.DataFrame(np.zeros(shape=(len(df0)*4,len(nams))),columns=nams)
    cnt=0
    for i in range(len(df0)):
        for j in range(4):
            tmp=np.isnan(df0['SPECIES_PCT_' + str(j+1)][i])
            if tmp==True:
                continue
            df1.iloc[cnt,0]=id[i]                              # Species rcrd id
            df1.iloc[cnt,1]=id[i]                              # layer rcrd id
            df1.iloc[cnt,2]=j+1                                # species id
            df1.iloc[cnt,3]=df0['SPECIES_CD_' + str(j+1)][i]   # species cd
            df1.iloc[cnt,4]=df0['SPECIES_PCT_' + str(j+1)][i]  # species pct
            df1.iloc[cnt,5]=df0['PROJ_AGE_1'][i]               # Age
            df1.iloc[cnt,6]=''                                 # Height %h_sps; % (include in BAH)
            df1.iloc[cnt,7]=df0['SITE_INDEX'][i]               # SI
            df1.iloc[cnt,8]=10                                 # ytbh
            df1.iloc[cnt,9]=12                                 # bhage
            df1.iloc[cnt,10]=''                                # site curve
            cnt=cnt+1
    df1.to_csv(pthout + '\\SPECIES.csv',sep=',',index=False)

    # Non-veg
    nams=['POLYGON_RCRD_ID','NON_VEG_ID','NON_VEG_COVER_TYPE','NON_VEG_COVER_PCT']
    df1=pd.DataFrame(np.zeros(shape=(len(df0),len(nams))),columns=nams)
    df1.iloc[:,0]=id                                # polygon id
    df1.iloc[:,1]=id                                # non veg id
    df1.iloc[:,2]='BR'                              # non veg type
    df1.iloc[:,3]=0                                 # non veg percent
    df1.to_csv(pthout + '\\NON_VEG.csv',sep=',',index=False)

    # History
    nams=['HISTORY_RCRD_ID','LAYER_RCRD_ID','HISTORY_ID','SILV_BASE_CD','START_YEAR','END_YEAR','LAYER_PCT']
    df1=pd.DataFrame(np.zeros(shape=(len(df0),len(nams))),columns=nams)
    df1.iloc[:,0]=id                                # History rcrd id
    df1.iloc[:,1]=id                                # Layer rcrd id
    df1.iloc[:,2]=id                                # History id
    df1.iloc[:,3]=''                                # SILV Base
    df1.iloc[:,4]=''                                # Start year ?
    df1.iloc[:,5]=''                                # End year ?
    df1.iloc[:,6]=''                                # Layer percent ?
    df1.to_csv(pthout + '\\HISTORY.csv',sep=',',index=False)

    # VRI Adjust
    nams=['VRIADJST_RCRD_ID','POLYGON_RCRD_ID','LAYER_ID','LH_075','BA_125','WSV_075','WSV_125','VCU_125','VD_125','VDW_125']
    df1=pd.DataFrame(np.zeros(shape=(len(df0),len(nams))),columns=nams)
    df1.iloc[:,0]=''                 # vriadj rcrd id
    df1.iloc[:,1]=''                 # polygon rcrd id
    df1.iloc[:,2]=1                                 # layer id
    df1.iloc[:,3]=''                                #
    df1.iloc[:,4]=''                                #
    df1.iloc[:,5]=''                                #
    df1.iloc[:,6]=''                                #
    df1.iloc[:,7]=''                                #
    df1.iloc[:,8]=''                                #
    df1.iloc[:,9]=''                                #
    df1.to_csv(pthout + '\\VRIADJST.csv',sep=',',index=False)

    return

#%%

def ProcessOutputsFromVDYP():

    # Input path to VDYP output data
    pthin=r'C:\Users\rhember\Documents\Data\VDYP7\FCI_DST_RehabFire\Outputs\tables.dat'

    # Output path
    pthout=r'C:\Users\rhember\Documents\Data\VDYP7\FCI_DST_RehabFire\Outputs'

    # Import data
    with open (pthin, "r") as myfile:
        data=myfile.readlines()


    # Initialize dataframe
    nams=['ID_Table','Year','Age','SPECIES_CD_1','SPECIES_PCT_1','Vws']
    df=pd.DataFrame(np.zeros(shape=(len(data),len(nams))),columns=nams)

    # Initialize table counter
    c_tab=0

    # Initiative counter
    cnt=0

    # Loop through tables
    for i in range(len(data)):

        # Update table counter
        if data[i][0:4]=='----':
            c_tab=c_tab+1

        if data[i][0]=='1' or data[i][0]=='2':

            # ID
            df['ID_Table'][cnt]=c_tab

            # Year
            df['Year'][cnt]=np.array(data[i][0:4]).astype(float)

            # Age
            df['Age'][cnt]=np.array(data[i][6:9]).astype(float)

            # Vws
            tmp=data[i][121:127]
            if tmp!='      ':
                df['Vws'][cnt]=np.array(tmp).astype(float)

            cnt=cnt+1


    # Save to file
    writer=pd.ExcelWriter(pthout + '\\VDYP_output_ForCBM3.xlsx')
    df.to_excel(writer,'Sheet1')
    writer.save()


    # Check to see if it worked
    ivl=np.arange(0,len(df),50)
    for i in range(6):
        ind=pl.find(df['ID_Table']==ivl[i])
        plt.plot(df['Age'][ind],df['Vws'][ind],'k.')
        #time.sleep(2)

    return

# Formatting of VDYP7 output:
#0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001111111111111111111111111111111111111111111111111111111111111
#0000000001111111111222222222223333333333444444444455555555556666666666777777777788888888889999999999000000000111111111122222222223333333333444444444455555555556
#1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
#Year  Age                      Stand Composition                      % Stk   SI   D Hgt  L Hgt   Dia    TPH       BA      Vws    Vcu    Vd     Vdw   Vdwb  Mode
#1908   50 FD  100.0       0.0       0.0       0.0       0.0       0.0  73.2  10.60   9.04   6.76  11.5   725.23   7.5278   18.7    7.1    7.1    7.1    7.0 Back




s