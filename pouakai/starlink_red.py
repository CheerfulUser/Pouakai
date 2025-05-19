from core import pouakai

from glob import glob
import numpy as np

path = '/home/phys/astro8/MJArchive/MOA/ALERT/'

files = glob(path + '*BW3*.gz')

for i in range(len(files)):
    print(files[i])
    #try:
    pouakai(file=files[i],savepath='/home/phys/astronomy/rri38/reduction/starlink/',time_tolerence=1000,dark_tolerence=60)
    #except:
        #print('failed ',files[i])
