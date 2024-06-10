from pouakai import consume_moa

from glob import glob
import numpy as np

path = '/home/phys/astro8/MJArchive/MOA/ALERT/2021/'
savepath = '/home/phys/astronomy/rri38/reduction/lsst_test/'

files = glob(path + '*lsst*.gz')
files = np.array(files)

print(files)
consume_moa(files,cores=10,savepath=savepath,rescale=True,update_cals=False)