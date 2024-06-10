from pouakai import consume_moa

from glob import glob
import numpy as np
import pandas as pd

path = '/home/phys/astro8/MJArchive/octans/'
savepath = '/home/phys/astronomy/rri38/reduction/BC/'

d = pd.read_csv('/home/phys/astronomy/rri38/fli_pouakai/pouakai/cal_lists/obs_list.csv')
#sn = '2022nah'
#sn = '2022rjg'
ind = []
for i in range(len(d)):
    file = d.iloc[i].filename
    if '2022nah'  in file:
        ind += [i]
    if '2022rjg' in file:
        ind += [i]

files = d.filename.values[ind]

print(files)
consume_moa(files,cores=5,savepath=savepath,rescale=True,update_cals=False,time_tolerence=200,dark_tolerence=0,plot=False)