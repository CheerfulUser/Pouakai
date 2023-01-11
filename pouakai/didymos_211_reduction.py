from core import pouakai

from glob import glob
import numpy as np

path = '/home/phys/astro8/MJArchive/MOA/ALERT/*/'

files = glob(path + '*idymos*3*.gz')
files = np.array(files)

for i in range(len(files)):
    print(files[i])
    try:
        pouakai(file=files[i],savepath='/home/phys/astronomy/rri38/reduction/dart/didymos/')
    except:
        print('failed ',files[i])
