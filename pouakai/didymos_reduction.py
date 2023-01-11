from core import pouakai
from joblib import Parallel, delayed
from glob import glob

path = '/home/phys/astro8/MJArchive/MOA/ALERT/*/'

files = glob(path + '*idymos*-3.fit.gz')
savepath='/home/phys/astronomy/rri38/moa/data/reduced/dart/didymos_no_zp/'
#pouakai(file=files[0],savepath=savepath)

def hack(file,savepath):
    #try:
    pouakai(file,savepath=savepath,rescale=False)
    #except Exception as e:
     #   print(e)

Parallel(n_jobs=30)(delayed(hack)(file,savepath) for file in files)
'''
for i in range(len(files)):    

    print(files[i])
    try:
        pouakai(file=files[i],savepath='/home/phys/astronomy/rri38/moa/data/reduced/dart/didymos/')
    except:
        print('failed ',files[i])
'''