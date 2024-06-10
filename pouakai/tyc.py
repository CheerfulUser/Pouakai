from pouakai import consume_moa

from glob import glob
import numpy as np
import pandas as pd
import time

start = time.time()

files = np.genfromtxt('/home/users/zgl12/tyc_files.csv', delimiter = ',', dtype = str)

savepath = '/home/users/zgl12/Test_Save_Pouakai/TYC/'

consume_moa(files,cores=1,savepath=savepath,rescale=True,update_cals=False,time_tolerence=20,dark_tolerence=0,plot=True, center = [112.021, -54.681])
print()
print(f"This took {(time.time()-start)/3600:.2f} hours")
