from pouakai import consume_moa
from glob import glob 
import numpy as np

files = glob('/home/phys/astro8/MJArchive/octans/**/EPS*')
files = np.array(files)
consume_moa(files=files,savepath='/home/phys/astronomy/rri38/fli/data/reduced/',cores=1,update_cals=False,overwrite=False,time_tolerence=1000)
#consume_moa(files=files,savepath='/home/phys/astronomy/rri38/fli/data/reduced/',cores=1,update_cals=False,overwrite=True)
