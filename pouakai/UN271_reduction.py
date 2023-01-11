from core import pouakai
from joblib import Parallel, delayed
from glob import glob
import numpy as np

path = '/home/phys/astro8/MJArchive/MOA/ALERT/*/'

files = glob(path + '*UN271*-3*.gz')
files = np.array(files)

savepath = '/home/phys/astronomy/rri38/moa/data/reduced/look/un271/'
def hack(file,savepath):
    #try:
    pouakai(file,savepath=savepath,rescale=True)
    #except Exception as e:
     #   print(e)

Parallel(n_jobs=2)(delayed(hack)(file,savepath) for file in files)

'''for i in range(len(files)):
    print(files[i])
    try:
        pouakai(file=files[i],savepath=)
    except:
        print('failed ',files[i])
'''