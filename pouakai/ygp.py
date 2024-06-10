from pouakai import consume_moa

from glob import glob
import numpy as np
import pandas as pd
import time

start = time.time()

files = np.genfromtxt('/home/users/zgl12/ygp_files.csv', delimiter = ',', dtype = str)

# files = files[13:]
# print(files[0])

# files = ['/home/phys/astro8/MJArchive/octans/20231123/2023ygp-0003_i.fit']
# print(files)

savepath = '/home/users/zgl12/Test_Save_Pouakai/YGP/'

consume_moa(files,cores=1,savepath=savepath,rescale=True,update_cals=False,time_tolerence=20,dark_tolerence=0,plot=False, center = [74.188, -68.555], limit_source=[12,17])
print()
print(f"This took {(time.time()-start)/3600:.2f} hours")
