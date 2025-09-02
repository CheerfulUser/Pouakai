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
package_directory = os.path.dirname(os.path.abspath(__file__)) + '/'
tmp = os.environ['TMPDIR']

class consume_moa():
    def __init__(self,files,savepath,time_tolerance=30,dark_tolerance=3,
				 local_astrom=True,verbose=True,rescale=True,update_cals=True,
                 cores=10, overwrite=False,compress=True,calibrate=True,
                 plot=True, limit_source = None, center = None, telescope = 'moa'):
        
        self.files = list(files)
        self.savepath = savepath
        self.verbose = verbose
        self._overwrite(overwrite=overwrite)
        self._clip_files()
        self.dark_tolerance = dark_tolerance
        self.time_tolerance = time_tolerance
        self.calibrate = calibrate
        self.compress = compress
        self.plot = plot
        self.center = center
        self.telescope = telescope
        self.limit_source = limit_source
        
        self.local_astrom = local_astrom
        self.rescale = rescale
        
        self.cores = cores
        self._kill_wcs_tmp()
        #running
        if update_cals:
            self._update_cals()
            self._update_masters()

        self.digest()
        self._kill_wcs_tmp()
        
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
        print('Updating cals')
        sort_cals(self.verbose, num_cores = self.cores, telescope = self.telescope)
    
    def _update_masters(self):
        print('Updating masters')
        make_masters(time_frame_dark = self.dark_tolerance, time_frame_flat = self.time_tolerance, num_cores=self.cores,verbose=self.verbose, telescope = self.telescope, dark_tolerance = 1)


    def _run_func(self,file):
        # process = psutil.Process()
        # print('start mem: ', process.memory_info().rss/1024**2)  # in bytes 
        
        p = pouakai(file,time_tolerance=self.time_tolerance,
                        dark_tolerance=self.dark_tolerance, savepath = self.savepath,
                        local_astrom=self.local_astrom,rescale=self.rescale,
                        verbose=self.verbose, center = self.center, limit_source = self.limit_source,
                        calibrate=self.calibrate,plot=self.plot, telescope = self.telescope)
        return p.log

    def _overwrite(self,overwrite):
        if not overwrite:
            #print('!!!!!!!!!!!!!!!!!!!!')
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
        ### ZAAAACCCC
        try:
            if self.telescope.lower() == 'moa':
                self.log = pd.read_csv(package_directory + 'cal_lists/moa_calibrated_image_list.csv')
            else:
                self.log = pd.read_csv(package_directory + 'cal_lists/bc_calibrated_image_list.csv')
        except:
            document = {'name':None, 'band':None,'chip':None, 'telescope':None, 'readout':None,
					'exptime':None,'jd':None,'date':None,
					'field':None,'filename':None,'flat':None,
					'dark':None,'tdiff_flat':None,
					'tdiff_dark':None,'zp':None,'zperr':None,
					'maglim5':None,'maglim3':None,'savename':None}

            self.log = pd.DataFrame(columns= document.keys())

    def _load_error_log(self):
        ##### ZAAAAACCCCC
        if self.telescope.lower() == 'moa':
            self.error_log = pd.read_csv(package_directory + 'cal_lists/moa_error_log.csv')
        else:
            self.error_log = pd.read_csv(package_directory + 'cal_lists/bc_error_log.csv')
    
    def _log_error(self,error):
        self._load_error_log()
        error 

    def _compress(self):
        compress = 'gzip -f ' + self.savepath+'/*.fits'
        os.system(compress)

    def _update_log(self):
        self._load_calibration_log()
        logs = glob(self.savepath + 'log/*.csv')
        # logs = glob(self.savepath + 'log/*.csv')
        for log in logs:
            new_entry = pd.read_csv(log)
            self.log = pd.concat([self.log, new_entry], ignore_index=True)
        
        if self.telescope.lower() == 'moa':
            self.log.to_csv(package_directory + 'cal_lists/moa_calibrated_image_list.csv',index=False)
        else:
            self.log.to_csv(package_directory + 'cal_lists/bc_calibrated_image_list.csv',index=False)
        os.system(f'rm -rf {self.savepath}log/*.csv')

    def digest(self):
        if (self.cores > 1) & (len(self.files) > 1):
            Parallel(self.cores)(delayed(self._run_func)(file) for file in self.files)
        else:
            for i in range(len(self.files)):
                self._run_func(self.files[i])

        self._update_log()

        #if self.compress:
         #   self._compress()