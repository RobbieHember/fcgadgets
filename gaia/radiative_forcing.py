'''
EXPLORING CALCULATION OF RADIATIVE FORCING FROM GHG EMISSIONS OUTPUT
'''

#%% Radiative forcing

gp=gu.SetGraphics('Manuscript')

# Import
v0=cbu.LoadSingleOutputFile(meta,1,0,0)
tv=np.arange(meta[pNam]['Project']['Year Start Saving'],meta[pNam]['Project']['Year End']+1,1)

# Extract biophysical parameters
bB=meta['Param']['BE']['Biophysical']

# Import atmospheric GHG concentration data
dA=gu.ReadExcel(meta[pNam]['Paths']['Model Code'] + '\\Parameters\\Parameters_RCP_GHG_Abundance_AR5.xlsx')
for k in dA.keys():
	if k=='Year':
		continue
	dA[k]=np.interp(tv,dA['Year'],dA[k])

#%% Radiative forcing

c_in=np.sum(v0['Atm_CO2_In'],axis=1)
c_out=np.sum(v0['Atm_CO2_Out'],axis=1)

t=np.arange(0,c_in.size)
r=0.001

c_out_decay=np.zeros( (c_in.size,c_in.size) )
for i in range(c_in.size):
    t=np.arange(0,c_in.size-i)
    c_out_decay[i:,i]=-np.append(0,np.diff(c_in[i]*np.exp(-r*t)))
c_out_decay=np.sum(c_out_decay,axis=1)

plt.close('all')
plt.plot(tv,c_in-c_out,'b-')
plt.plot(tv,c_in-c_out-c_out_decay,'r--')

#%% Radiative forcing

for i in range(c_in.size):
    #plt.plot(np.diff(c_out_decay))
    c_atm[i:,i]=c_atm[i:,i]+c_out_decay

# Remove outputs
#c_atm=np.mean(c_atm,axis=1)#-c_out

plt.close('all')
plt.plot(tv,c_atm,'b-')
plt.plot(tv,np.cumsum(c_in-c_out),'r--')
#plt.plot(tv,c_atm-c_out,'r--')

#plt.close('all')
#plt.plot(tv,c_in-c_out,'r-')
#plt.plot(tv,np.cumsum(c_out),'g-')

#%% Radiative forcing
#
plt.plot(c1[:,0]/c1[0,0],'b-')
#plt.plot(a[:,0],'b-')


#bB['Ratio_C_to_CO2']*
#/bB['Atmospheric Density CO2']