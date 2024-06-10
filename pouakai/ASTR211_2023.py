from pouakai import consume_moa

from glob import glob
import numpy as np
import pandas as pd

path = '/home/phys/astro8/MJArchive/octans/'
savepath = '/home/phys/astronomy/rri38/reduction/BC/'

d = pd.read_csv('/home/phys/astronomy/rri38/fli_pouakai/pouakai/cal_lists/obs_list.csv')
dates = d['date'].values

ind = []

for i in range(len(dates)):
    if ('2023-08' in dates[i]) | ('2023-09' in dates[i]):
        ind += [i]

dd = d.iloc[ind]

dd = dd[~dd['band'].str.contains('bad')]

ind = []

for i in range(len(dd)):
    if 'gaia712' in dd['object'].values[i].lower():
        ind += [i]
dd = dd.iloc[ind]

files = dd.filename.values

print(files)
consume_moa(files,cores=5,savepath=savepath,rescale=True,update_cals=False,time_tolerence=200,dark_tolerence=0,plot=False)