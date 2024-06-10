from pouakai import consume_moa

from glob import glob
import numpy as np
import pandas as pd

path = '/home/phys/astro8/MJArchive/octans/'
savepath = '/home/phys/astronomy/rri38/reduction/BC/'

d = [20231117,20231118,20231120,20231121]
#sn = '2022nah'
#sn = '2022rjg'
runner = []
for i in range(len(d)):
    files = glob(path+str(d[i])+'/*')
    for file in files:
        if 'C2021'  in file:
            runner += [file]

print(runner)

print(runner)
consume_moa(runner,cores=5,savepath=savepath,rescale=True,update_cals=False,time_tolerence=200,dark_tolerence=0,plot=False)