from pouakai import consume_moa
from glob import glob
import numpy as np

savepath = '/home/phys/astronomy/rri38/reduction/bw3/'
path = '/home/phys/astro8/MJArchive/MOA/ALERT/'

files = glob(path + '*WB3*.fit.gz')
files = np.array(files)

consume_moa(files,cores=10,savepath=savepath,rescale=True,update_cals=False,time_tolerence=300,overwrite=False)
