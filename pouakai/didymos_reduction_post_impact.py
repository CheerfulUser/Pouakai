from pouakai import consume_moa
#from joblib import Parallel, delayed
from glob import glob
import numpy as np

path = '/home/phys/astro8/MJArchive/MOA/ALERT/'

files = glob(path + '*idymos*-R-3.fit.gz')#np.append(glob(path + '2022/*idymos*-R-3.fit.gz'),glob(path + '*idymos*-R-3.fit.gz'))
savepath='/home/phys/astronomy/rri38/moa/data/reduced/dart/didymos/night1/no_cor/'

anum = []
for file in files:
    anum += [file.split('/')[-1].split('-')[0]]
done = glob(savepath + 'red/*.gz')
dnum =  []
for file in done:
    dnum += [file.split('/')[-1].split('-')[0]]

not_done = set(dnum) ^ set(anum)
not_done = list(not_done)
not_done.sort()

todo = []
for d in not_done:
    todo += glob(path + d + '-*idymos*-R-3.fit.gz')

print('!!! to do: ',len(todo))
todo.sort()

consume_moa(todo,cores=60,savepath=savepath,rescale=False)

#pouakai(file=files[0],savepath=savepath)
'''
def hack(file,savepath):
    #try:
    pouakai(file,savepath=savepath,rescale=False)
    #except Exception as e:
     #   print(e)

Parallel(n_jobs=20)(delayed(hack)(file,savepath) for file in todo)
'''
print('!!!DONE!!!')
'''
for i in range(len(files)):    

    print(files[i])
    try:
        pouakai(file=files[i],savepath='/home/phys/astronomy/rri38/moa/data/reduced/dart/didymos/')
    except:
        print('failed ',files[i])
'''