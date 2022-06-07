
#------------------------------------------------------------------------------
# Read geodatabase
#------------------------------------------------------------------------------

import geopandas as gpd
import fiona
from bokeh.plotting import figure
from bokeh.io import output_notebook, show

fin=r"E:\Data\ForestInventory\PSP-NADB\Data\01_RawFiles\QC\Release_2017-07\Pep\PEP.gdb"

fiona.listlayers(fin)

dat=gpd.read_file(fin, layer='PLACETTE_MES')
sorted(list(dat.columns.values))
dat.to_csv(r"E:\Data\ForestInventory\PSP-NADB\Data\01_RawFiles\QC\Release_2017-07\PLACETTE_MES.csv")

dat=gpd.read_file(fin, layer='PLACETTE')
sorted(list(dat.columns.values))
dat.to_csv(r"E:\Data\ForestInventory\PSP-NADB\Data\01_RawFiles\QC\Release_2017-07\PLACETTE.csv")

dat=gpd.read_file(fin, layer='STATION_PE')
sorted(list(dat.columns.values))
dat.to_csv(r"E:\Data\ForestInventory\PSP-NADB\Data\01_RawFiles\QC\Release_2017-07\STATION_PE.csv")

dat=gpd.read_file(fin, layer='DENDRO_ARBRES')
sorted(list(dat.columns.values))
dat.to_csv(r"E:\Data\ForestInventory\PSP-NADB\Data\01_RawFiles\QC\Release_2017-07\DENDRO_ARBRES.csv") # , header=False
L=list(dat.columns.values)
L.to_csv(r"E:\Data\ForestInventory\PSP-NADB\Data\01_RawFiles\QC\Release_2017-07\DENDRO_ARBRES_lab.csv")

dat=gpd.read_file(fin, layer='DENDRO_ARBRES_ETUDES')
sorted(list(dat.columns.values))
dat.to_csv(r"E:\Data\ForestInventory\PSP-NADB\Data\01_RawFiles\QC\Release_2017-07\DENDRO_ARBRES_ETUDES.csv")

