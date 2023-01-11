from core import pouakai
from sort_images import sort_cals
from calibration_masters import make_masters
from joblib import Parallel, delayed
import pandas as pd
import numpy as np

import os
package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'

class consume_moa():
    def __init__(self,files,savepath,time_tolerence=100,dark_tolerence=10,
				 local_astrom=True,verbose=True,rescale=True,update_cals=True,cores=10):
        
        self.files = list(files)
        self._clip_files()
        self.dark_tolerence = dark_tolerence
        self.time_tolerence = time_tolerence
        self.savepath = savepath
        self.local_astrom = local_astrom
        self.rescale = rescale
        self.verbose = verbose
        self.cores = cores

        #running
        if update_cals:
            self._update_cals()
            self._update_masters()

        self.digest()
        
    def _clip_files(self):
        ind = []
        for i in range(len(self.files)):
            if '-99.fit.gz' not in self.files[i]:
                if '-0.fit.gz' not in self.files[i]:
                    ind += [i]

        self.files = list(np.array(self.files)[ind])
    
    def _update_cals(self):
        sort_cals(self.verbose)
    
    def _update_masters(self):
        make_masters(self.verbose)


    def _run_func(self,file):
        try:
            p = pouakai(file,time_tolerence=self.time_tolerence,
                            dark_tolerence=self.dark_tolerence, savepath = self.savepath,
                            local_astrom=self.local_astrom,rescale=self.rescale,verbose=self.verbose)
            return p.log
        except Exception as e:
            
            #self._log_error(e)
            print('Failed')


    def _load_calibration_log(self):
        self.log = pd.read_csv(package_directory + 'cal_lists/calibrated_image_list.csv')

    def _load_error_log(self):
        self.error_log = pd.read_csv(package_directory + 'cal_lists/error_log.csv')
    
    def _log_error(self,error):
        self._load_error_log()
        error 


    def digest(self):
        if (self.cores > 1) & (len(self.files) > 1):
            entries = Parallel(self.cores)(delayed(self._run_func)(file) for file in self.files)
        else:
            entries = []
            for i in range(len(self.files)):
                entries += [self._run_func(self.files[i])]

        self._load_calibration_log()

        for entry in entries:
            new_entry = pd.DataFrame([entry])
            self.log = pd.concat([self.log, new_entry], ignore_index=True)
		
        self.log.to_csv(package_directory + 'cal_lists/calibrated_image_list.csv',index=False)

        self._load_calibration_log()