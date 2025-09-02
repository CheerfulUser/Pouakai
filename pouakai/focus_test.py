from pouakai import consume_moa
from glob import glob
import numpy as np


path = '/home/phys/astro8/MJArchive/MOA/ALERT/'

files = glob(path + '*focus_NGC2451*.gz')
files = np.array(files)

ind = []
for i in range(len(files)):
    if '-99.fit.gz' not in files[i]:
        if '-0.fit.gz' not in files[i]:
            ind += [i]
files = files[ind]

print('!!!! ',len(files))
savepath = '/home/phys/astronomy/rri38/moa/data/reduced/focus_test/'
consume_moa(files,cores=2,savepath=savepath,rescale=False,update_cals=False)