from pouakai import consume_moa

from glob import glob
import numpy as np

path = '/home/phys/astro8/MJArchive/MOA/ALERT/'
savepath = '/home/phys/astronomy/rri38/reduction/astr211/2023/'

f = glob(path+'A1[5-6][0-9][0-9][0-9]*-3.fit.gz')
rn = [x.split('/')[-1] for x in f]
rn = [int(x.split('-')[0].split('A')[1]) for x in rn]
rn = np.array(rn)
ind = rn > 15264
rn = rn[ind]
files = np.array(f)[ind]

#red = glob('/home/phys/astronomy/rri38/reduction/astr211/2023/*.gz')
#redrn = [x.split('/')[-1] for x in red]
#redrn = [int(x.split('-')[0].split('A')[1]) for x in redrn]

#indo = []
#for i in range(len(rn)):
#    if rn[i] not in redrn:
#        indo += [i]

#files = files[indo]

print(files)
consume_moa(files,cores=10,savepath=savepath,rescale=True,update_cals=False,time_tolerence=1e5,dark_tolerence=60,overwrite=False,plot=False)