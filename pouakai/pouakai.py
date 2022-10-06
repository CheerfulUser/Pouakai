from core import pouakai
from sort_images import sort_cals
from calibration_masters import make_masters
from joblib import Parallel, delayed
import pandas as pd

import os
package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'

class comsume_moa():
    def __init__(self,files,time_tolerence=100,dark_tolerence=10,savepath='',
				 local_astrom=True,verbose=True,rescale=True,update_cals=True,cores=10):
        
        self.files = list(files)
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
        

    
    def _update_cals(self):
        sort_cals(self.verbose)
    
    def _update_masters(self):
        make_masters(self.verbose)


    def _run_func(self,file):
        log = pouakai(file,time_tolerence=self.time_tolerence,
                        dark_tolerence=self.dark_tolerence, savepath = self.savepath,
                        local_astrom=self.local_astrom,rescale=self.rescale,verbose=self.verbose,
                        output_record=True)
        return log


    def _load_calibration_log(self):
        self.log = pd.read_csv(package_directory + 'cal_lists/calibrated_image_list.csv')

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