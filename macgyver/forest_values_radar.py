'''
FOREST VALUES
'''

#%% Modules
import matplotlib.pyplot as plt
import numpy as np
import fcgadgets.macgyver.util_general as gu

#%% Business-as-usual

labels=['Timber\nProduction','Climate\nMitigation','Forest\nConservation','Air\nQuality','Water\nSupply','Grassland\nConservation','Wildlife\nHabitat','Economic\nDiversity', \
		'Terrain\nStability','Pasture\nHealth','Water\nQuality','Tourism','Wildfire\nMitigation','Riparian &\nWetland\nConservation','Flood\nMitigation']

values=[5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]
num_vars=len(labels)

# Split the circle into even parts and save the angles
# so we know where to put each axis.
angles=np.linspace(0,2*np.pi, num_vars, endpoint=False).tolist()

#values += values[:1]
#angles += angles[:1]
values.insert(0,values[0])
angles.insert(0,angles[0])
labels.insert(0,labels[0])

plt.close('all');fig,ax=plt.subplots(figsize=(gu.cm2inch(12,12)),subplot_kw=dict(polar=True)); fn=10
ax.plot(angles,values,'o-',linewidth=0.5,ms=2,color=[0.24,0.49,0.79])
ax.fill(angles,values,alpha=0.25,color=[0.24,0.49,0.79])
ax.set_theta_offset(np.pi/2) # Fix axis to go in the right order and start at 12 o'clock.
ax.set_theta_direction(-1)
ax.set_thetagrids(np.degrees(angles),labels,fontsize=fn,fontweight='bold') # Draw axis lines for each angle and label.
#ax.grid(True)
ax.set_rlim(0,10)
ax.set_rlabel_position(62)
plt.yticks(np.arange(1,10,1),color="grey",size=9)
plt.tight_layout()
#ax.spines['inner'].set_color('w')

path=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Mitigation'
gu.PrintFig(path + '\\ForestValues_Landscape1','png',900)

#%% Underplanting

#labels=['Timber\nProduction','Climate\nMitigation','Forest\nConservation','Air\nQuality','Water\nSupply','Grassland\nConservation','Wildlife\nHabitat','Economic\nDiversity', \
#		'Terrain\nStability','Pasture\nHealth','Water\nQuality','Tourism','Wildfire\nMitigation','Riparian &\nWetland\nConservation','Flood\nMitigation']
values2=[7,8,5,5,6,5,5,6,8,5,6,5,4,5,7]

#values2 += values2[:1]
values2.insert(0,values2[0])

plt.close('all');fig,ax=plt.subplots(figsize=(gu.cm2inch(12,12)),subplot_kw=dict(polar=True))
ax.plot(angles,values2,'o-',linewidth=0.5,ms=2,color=[0.6,0.9,0])
ax.fill(angles,values2,alpha=0.25,color=[0.6,0.9,0])
ax.plot(angles,values,'-',linewidth=0.5,ms=2,color=[0.24,0.49,0.79])
ax.fill(angles,values,alpha=0.25,color=[0.24,0.49,0.79])
ax.set_theta_offset(np.pi/2) # Fix axis to go in the right order and start at 12 o'clock.
ax.set_theta_direction(-1)
ax.set_thetagrids(np.degrees(angles),labels,fontsize=fn,fontweight='bold') # Draw axis lines for each angle and label.
#ax.grid(True)
ax.set_rlim(0,10)
ax.set_rlabel_position(62)
plt.yticks(np.arange(1,10,1),color="grey",size=9)
plt.tight_layout()
#ax.spines['inner'].set_color('w')

path=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Mitigation'
gu.PrintFig(path + '\\ForestValues_Underplanting','png',900)

#%% Forest Nutrient Management

#labels=['Timber\nProduction','Climate\nMitigation','Forest\nConservation','Air\nQuality','Water\nSupply','Grassland\nConservation','Wildlife\nHabitat','Economic\nDiversity', \
#		'Terrain\nStability','Pasture\nHealth','Water\nQuality','Tourism','Wildfire\nMitigation','Riparian &\nWetland\nConservation','Flood\nMitigation']
values2=[7,8,5,5,5,5,5,5,5,5,5,5,5,5,5]

values2.insert(0,values2[0])

plt.close('all');fig,ax=plt.subplots(figsize=(gu.cm2inch(12,12)),subplot_kw=dict(polar=True))
ax.plot(angles,values2,'o-',linewidth=0.5,ms=2,color=[0.6,0.9,0])
ax.fill(angles,values2,alpha=0.25,color=[0.6,0.9,0])

ax.plot(angles,values,'-',linewidth=0.5,ms=2,color=[0.24,0.49,0.79])
ax.fill(angles,values,alpha=0.25,color=[0.24,0.49,0.79])

ax.set_theta_offset(np.pi/2) # Fix axis to go in the right order and start at 12 o'clock.
ax.set_theta_direction(-1)
ax.set_thetagrids(np.degrees(angles),labels,fontsize=fn,fontweight='bold') # Draw axis lines for each angle and label.
#ax.grid(True)
ax.set_rlim(0,10)
ax.set_rlabel_position(62)
plt.yticks(np.arange(1,10,1),color="grey",size=9)
plt.tight_layout()
#ax.spines['inner'].set_color('w')

path=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Mitigation'
gu.PrintFig(path + '\\ForestValues_FNM','png',900)

#%% Intensified management corridors

values2=[5,7,7,7,5,5,5,7,5,5,5,5,7,7,5]

values2.insert(0,values2[0])

plt.close('all');fig,ax=plt.subplots(figsize=(gu.cm2inch(12,12)),subplot_kw=dict(polar=True))
ax.plot(angles,values2,'o-',linewidth=0.5,ms=2,color=[0.6,0.9,0])
ax.fill(angles,values2,alpha=0.25,color=[0.6,0.9,0])

ax.plot(angles,values,'-',linewidth=0.5,ms=2,color=[0.24,0.49,0.79])
ax.fill(angles,values,alpha=0.25,color=[0.24,0.49,0.79])

#ax.set_theta_offset(np.pi/2) # Fix axis to go in the right order and start at 12 o'clock.
#ax.set_theta_direction(-1)
ax.set_thetagrids(np.degrees(angles),labels,fontsize=fn,fontweight='bold') # Draw axis lines for each angle and label.
#ax.grid(True)
ax.set_rlim(0,10)
ax.set_rlabel_position(62)
plt.yticks(np.arange(1,10,1),color="grey",size=9)
plt.tight_layout()
#ax.spines['inner'].set_color('w')

path=r'C:\Users\rhember\OneDrive - Government of BC\Figures\Mitigation'
gu.PrintFig(path + '\\ForestValues_Landscape2','png',900)

#%%
