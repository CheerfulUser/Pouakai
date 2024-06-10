from core import pouakai
from sort_images import sort_cals
from calibration_masters import make_masters
from joblib import Parallel, delayed
import pandas as pd
import numpy as np
from glob import glob
import os
import gc
from copy import deepcopy
package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'
#tmp = os.environ['TMPDIR']
tmp = '/home/users/zgl12/Temp_Dir/'

class consume_moa():
    def __init__(self,files,savepath,time_tolerence=60,dark_tolerence=1,
				 local_astrom=True,verbose=True,rescale=True,update_cals=True,
                 cores=1, overwrite=False,plot=False,center=None,limit_source = None):
        
        self.savepath = savepath
        self.files = list(files)
	#print(savepath)
        self.verbose = verbose
        self._overwrite(overwrite=overwrite)
        self._clip_files()
        self.dark_tolerence = dark_tolerence
        self.time_tolerence = time_tolerence
        self.limit_source = limit_source
        
        self.local_astrom = local_astrom
        self.rescale = rescale
        self.plot = plot
        self.center = center
        
        self.cores = cores

        #running
        if update_cals:
            self._update_cals()
            self._update_masters()
        self.digest()
        self._remove_tmp()

    def _remove_tmp(self):
        call = f'rm -rf {tmp}'
        os.system(call)
        
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
                            local_astrom=self.local_astrom,rescale=self.rescale,verbose=self.verbose,plot=self.plot,center=self.center,limit_source=self.limit_source)
            return p.log

        except Exception as e:
            self._log_error(e)
            print('!!! Failed: ',file)
            print(e)

    def _overwrite(self,overwrite):
        if not overwrite:
            try: 
                anum = []
                for file in self.files:
                    anum += [file.split('/')[-1].split('.')[0]]
                done = glob(self.savepath + 'cal/*.gz')
                dnum =  []
                for file in done:
                    dnum += [file.split('/')[-1].split('_')[0]]
                dnum = set(dnum)
                todo = []
                for i in range(len(anum)):
                    if not dnum.intersection({anum[i]}):
                        todo += [i]
                todo = np.array(todo)

                if self.verbose:
                    print(f'Dropping {len(anum) - len(todo)} files that have already been processed')
                self.files = list(np.array(self.files)[todo])
            except:
                pass



    def _load_calibration_log(self):
        self.log = pd.read_csv(package_directory + 'cal_lists/calibrated_image_list.csv')

    def _load_error_log(self):
        self.error_log = pd.read_csv(package_directory + 'cal_lists/error_log.csv')
    
    def _log_error(self,error):
        self._load_error_log()
        error 


    def digest(self):
        if (self.cores > 1) & (len(self.files) > 1):
            Parallel(self.cores)(delayed(self._run_func)(file) for file in self.files)
        else:
            for i in range(len(self.files)):
                self._run_func(self.files[i])
        self._update_log()
        
        

    def _update_log(self):
        self._load_calibration_log()
        logs = glob(self.savepath + 'log/*.csv')
        for log in logs:
            new_entry = pd.read_csv(log)
            self.log = pd.concat([self.log, new_entry], ignore_index=True)
        
        self.log.to_csv(package_directory + 'cal_lists/calibrated_image_list.csv',index=False)
        os.system(f'rm -rf {self.savepath}log/*.csv')
