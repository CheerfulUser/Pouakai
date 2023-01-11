from core import pouakai
from joblib import Parallel, delayed
from glob import glob
import numpy as np

def hack(file,savepath):
    #try:
    pouakai(file,savepath=savepath,rescale=True)
    #except Exception as e:
     #   print(e)

path = '/home/phys/astro8/MJArchive/MOA/ALERT/*/'

files = glob(path + '*C2017*-3*.gz')
files = np.array(files)
savepath = '/home/phys/astronomy/rri38/moa/data/reduced/look/c2017_w2/'

Parallel(n_jobs=10)(delayed(hack)(file,savepath) for file in files)


files = glob(path + '*C2019*-3*.gz')
files = np.array(files)
savepath = '/home/phys/astronomy/rri38/moa/data/reduced/look/c2019_e3/'
pouakai(files[0],savepath=savepath,rescale=True)
Parallel(n_jobs=10)(delayed(hack)(file,savepath) for file in files)

files = glob(path + '*C2021_S3*-3*.gz')
files = np.array(files)
savepath = '/home/phys/astronomy/rri38/moa/data/reduced/look/c2021_s3/'
pouakai(files[0],savepath=savepath,rescale=True)
Parallel(n_jobs=10)(delayed(hack)(file,savepath) for file in files)



'''for i in range(len(files)):
    print(files[i])
    try:
        pouakai(file=files[i],savepath=)
    except:
        print('failed ',files[i])
'''