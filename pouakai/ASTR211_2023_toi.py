from pouakai import consume_moa

from glob import glob
import numpy as np

path = '/home/phys/astro8/MJArchive/MOA/ALERT/2021/'
savepath = '/home/phys/astronomy/rri38/reduction/astr211/2023/'

files = glob(path+'*K2_45b*-3.fit.gz')
files = np.array(files)
rn = [int(x.split('/A')[-1].split('-')[0]) for x in files]
rn = np.array(rn)


red = glob('/home/phys/astronomy/rri38/reduction/astr211/2023/*.gz')
redrn = [x.split('/')[-1] for x in red]
redrn = [int(x.split('-')[0].split('A')[1]) for x in redrn]

indo = []
for i in range(len(rn)):
    if rn[i] not in redrn:
        indo += [i]

files = files[indo]

print(files)
consume_moa(files,cores=5,savepath=savepath,rescale=True,update_cals=False,time_tolerence=1e5,dark_tolerence=60,overwrite=False)