from pouakai import consume_moa
from glob import glob
import numpy as np

path = '/home/phys/astro8/MJArchive/MOA/ALERT/2021/'
savepath = '/home/phys/astronomy/rri38/reduction/dart_zp_test/'

path = '/home/phys/astro8/MJArchive/MOA/ALERT/'

files = np.append(glob(path + '2022/*idymos*.fit.gz'),glob(path + '*idymos*.fit.gz'))


print(files)
consume_moa(files,cores=1,savepath=savepath,rescale=True,update_cals=False,overwrite=False,verbose=False)
