from core import pouakai

from glob import glob

path = '/home/phys/astro8/MJArchive/MOA/ALERT/'

files = glob(path + '*12932*-3.fit.gz')

for i in range(len(files)):
    print(files[i])
    pouakai(file=files[i],savepath='/home/phys/astronomy/rri38/test/')