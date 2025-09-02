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

files = glob(path + '*SN2022rgz*-3*.gz')




f = glob(path + '*SN2022r_g_z*-3*.gz')
f = np.array(f)
files = np.append(files,f,axis=0)
f = glob(path + '*SN2022rjg*-3*.gz')
f = np.array(f)
files = np.append(files,f,axis=0)
f = glob(path + '*SN2022rng*-3*.gz')
f = np.array(f)
files = np.append(files,f,axis=0)
f = glob(path + '*SN2022rwq*-3*.gz')
f = np.array(f)
files = np.append(files,f,axis=0)
f = glob(path + '*SN2022rpl*-3*.gz')
f = np.array(f)
files = np.append(files,f,axis=0)

print('!!!! ',len(files))
savepath = '/home/phys/astronomy/rri38/moa/data/reduced/211SN/'
Parallel(n_jobs=20)(delayed(hack)(file,savepath) for file in files)



'''for i in range(len(files)):
    print(files[i])
    try:
        pouakai(file=files[i],savepath=)
    except:
        print('failed ',files[i])
'''