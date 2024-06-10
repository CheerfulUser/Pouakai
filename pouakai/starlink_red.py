from core import pouakai

from glob import glob
import numpy as np

path = '/home/phys/astro8/MJArchive/MOA/ALERT/'

files = glob(path + '*starlink*.fit.gz')
files = np.array(files)

ind = []
for i in range(len(files)):
    if '-99.fit.gz' not in files[i]:
        if '-0.fit.gz' not in files[i]:
            ind += [i]

files = files[ind]

for i in range(len(files)):
    print(files[i])
    #try:
    pouakai(file=files[i],savepath='/home/phys/astronomy/rri38/reduction/starlink/',time_tolerence=1000,dark_tolerence=60)
    #except:
        #print('failed ',files[i])
