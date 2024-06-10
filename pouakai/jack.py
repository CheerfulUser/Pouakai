from pouakai import consume_moa

import os
from glob import glob
import numpy as np
import pandas as pd
import time

start = time.time()

files = np.genfromtxt('/home/users/zgl12/jack_files.csv', delimiter = ',', dtype = str)

# files = files[13:]
# print(files[0])

# files = ['/home/phys/astro8/MJArchive/octans/20231123/2023ygp-0003_i.fit']
# print(files)

savepath = '/home/users/zgl12/Image_Processing/Test_Save_Pouakai/Jack/'

if not os.path.exists(savepath):
    os.makedirs(savepath)

consume_moa(files,cores=1,savepath=savepath,rescale=True,update_cals=False,time_tolerence=45,dark_tolerence=1,plot=False, limit_source=[12,17])
print()
print(f"This took {(time.time()-start)/3600:.2f} hours")