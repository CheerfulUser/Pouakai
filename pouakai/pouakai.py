from core import pouakai
from sort_images import sort_cals
from calibration_masters import make_masters
from joblib import Parallel, delayed
import pandas as pd
import numpy as np
from glob import glob
import os, psutil
import gc
package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'

class consume_moa():
    def __init__(self,files,savepath,time_tolerence=100,dark_tolerence=10,
				 local_astrom=True,verbose=True,rescale=True,update_cals=True,
                 cores=10, overwrite=False,compress=True,calibrate=True,plot=True):
        
        self.files = list(files)
        self.savepath = savepath
        self.verbose = verbose
        self._overwrite(overwrite=overwrite)
        self._clip_files()
        self.dark_tolerence = dark_tolerence
        self.time_tolerence = time_tolerence
        self.calibrate = calibrate
        self.compress = compress
        self.plot = plot
        
        self.local_astrom = local_astrom
        self.rescale = rescale
        
        self.cores = cores
        self._kill_wcs_tmp()
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

    def _kill_wcs_tmp(self):
        call = f'rm -rf {self.savepath}red/wcs_tmp/'
        os.system(call)

    def _update_cals(self):
        sort_cals(self.verbose)
    
    def _update_masters(self):
        make_masters(self.verbose)


    def _run_func(self,file):
        process = psutil.Process()
        print('start mem: ', process.memory_info().rss/1024**2)  # in bytes 
        
        try:
            p = pouakai(file,time_tolerence=self.time_tolerence,
                            dark_tolerence=self.dark_tolerence, savepath = self.savepath,
                            local_astrom=self.local_astrom,rescale=self.rescale,verbose=self.verbose,
                            calibrate=self.calibrate,plot=self.plot)
            from ctypes import cdll, CDLL
            cdll.LoadLibrary('libc.so.6')
            libc = CDLL('libc.so.6')
            libc.malloc_trim(0)
        
        except Exception as e:
         #   self._log_error(e)
            print('Failed: ' + file)
            print(e)

    def _overwrite(self,overwrite):
        if not overwrite:
            #print('!!!!!!!!!!!!!!!!!!!!')
            #try: 
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
                print(f'Droping {len(anum) - len(todo)} files that have already been processed')
            self.files = list(np.array(self.files)[todo])
            #except:
             #   pass



    def _load_calibration_log(self):
        self.log = pd.read_csv(package_directory + 'cal_lists/calibrated_image_list.csv')

    def _load_error_log(self):
        self.error_log = pd.read_csv(package_directory + 'cal_lists/error_log.csv')
    
    def _log_error(self,error):
        self._load_error_log()
        error 


    def _compress(self):
        compress = 'gzip -f ' + self.savepath+'/*.fits'
        os.system(compress)

    def _update_log(self):
        self._load_calibration_log()
        logs = glob(self.savepath + 'log/*.csv')
        for log in logs:
            new_entry = pd.read_csv(log)
            self.log = pd.concat([self.log, new_entry], ignore_index=True)
        
        self.log.to_csv(package_directory + 'cal_lists/calibrated_image_list.csv',index=False)
        os.system(f'rm -rf {self.savepath}log/*.csv')

    def digest(self):
        if (self.cores > 1) & (len(self.files) > 1):
            Parallel(self.cores)(delayed(self._run_func)(file) for file in self.files)
        else:
            for i in range(len(self.files)):
                self._run_func(self.files[i])

        self._load_calibration_log()

        self._update_log()

        #if self.compress:
         #   self._compress()