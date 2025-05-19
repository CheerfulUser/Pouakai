from pouakai import consume_moa

from glob import glob
import numpy as np

path = '/home/phys/astro8/MJArchive/MOA/ALERT/2022/'
savepath = '/home/phys/astronomy/rri38/reduction/sn/'

files = glob(path + '*22S05*-3.fit.gz') + glob(path + '*22S02*-3.fit.gz')
files = np.array(files)

print(files)
consume_moa(files,cores=1,savepath=savepath,rescale=True,update_cals=False,time_tolerence=1e5,dark_tolerence=60,overwrite=True)