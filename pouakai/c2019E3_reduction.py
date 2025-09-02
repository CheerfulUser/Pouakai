from pouakai import consume_moa
from glob import glob
import numpy as np

savepath = '/home/phys/astronomy/rri38/reduction/c2019E3/'
path = '/home/phys/astro8/MJArchive/MOA/ALERT/'

files = glob(path + '*c2019E3-R-3.fit.gz')
files = np.array(files)

consume_moa(files,cores=10,savepath=savepath,rescale=True,update_cals=False,time_tolerence=1000,overwrite=False,calibrate=False)
