from fooof import FOOOFGroup
import numpy as np
import scipy.io
import os
import time

v=33

SUBJLIST = [4,5,6,7,8,9,10,11,12,13,15,16,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34]

for isubj in SUBJLIST:
    for iblock in range(1,3):
        for m in range(1,4):

            if os.path.exists('/home/tpfeffer/pupmod/proc/src/pupmod_rest_fooof_s%d_m%d_b%d_v%d_proc.txt' % (isubj,m,iblock,v)) == True:
                continue

            os.system('touch /home/tpfeffer/pupmod/proc/src/pupmod_rest_fooof_s%d_m%d_b%d_v%d_proc.txt' % (isubj,m,iblock,v))
            try:
                tmp = scipy.io.loadmat('/home/tpfeffer/pupmod/proc/src/pupmod_rest_peakfreq_s%d_m%d_b%d_v%d.mat' % (isubj,m,iblock,v))
                dat = {}
                dat['fxx'] = tmp['outp'][0]['fxx'][0]
                dat['pxx'] = tmp['outp'][0]['pxx'][0]
            except:
                print("Error: File not found!")
                continue

            print('Processing S%d B%d M%d ...' % (isubj,iblock,m))


            freqs = np.squeeze(dat['fxx'])
            aper = np.empty([2,dat['pxx'].shape[1]])
            only_gauss = np.zeros([dat['pxx'].shape[0],dat['pxx'].shape[1]])
            full_gauss = np.zeros([dat['pxx'].shape[0],dat['pxx'].shape[1]])


            fm = FOOOFGroup(peak_width_limits=[1, 8], min_peak_height=0.05, max_n_peaks=6)                
            fm._maxfev = 30000
            freq_range = [3, 40]
            fm.fit(freqs, np.transpose(dat['pxx']), freq_range)
            tmp = fm.get_results()

            for isens in range(0,dat['pxx'].shape[1]):
              aper[:,isens] = tmp[isens].aperiodic_params

            F = fm.freqs

            for isens in range(0,dat['pxx'].shape[1]):
              for i in range(0,len(tmp[isens].gaussian_params)):
                  c = tmp[isens].gaussian_params[i][0]
                  w = tmp[isens].gaussian_params[i][2]
                  a = tmp[isens].gaussian_params[i][1]
                  only_gauss[:,isens] = only_gauss[:,isens] + a * np.exp ((-(F-c)**2)/(2*pow(w,2))) 
              b = tmp[isens].aperiodic_params[0]
              e = tmp[isens].aperiodic_params[1]
              full_gauss[:,isens] = only_gauss[:,isens] + b - np.log10(pow(F,e)) 

            scipy.io.savemat('/home/tpfeffer/pupmod/proc/src/pupmod_rest_fooof_s%d_m%d_b%d_v%d.mat' % (isubj,m,iblock,v), {'only_gauss': only_gauss, 'full_gauss': full_gauss,  'aper':  aper})

    
for isubj in SUBJLIST:
    for iblock in range(1,3):
        for m in range(1,4):

            if os.path.exists('/home/tpfeffer/pupmod/proc/src/pupmod_task_fooof_s%d_m%d_b%d_v%d_proc.txt' % (isubj,m,iblock,v)) == True:
                continue

            os.system('touch /home/tpfeffer/pupmod/proc/src/pupmod_task_fooof_s%d_m%d_b%d_v%d_proc.txt' % (isubj,m,iblock,v))
            try:
                tmp = scipy.io.loadmat('/home/tpfeffer/pupmod/proc/src/pupmod_task_peakfreq_s%d_m%d_b%d_v%d.mat' % (isubj,m,iblock,v))
                dat = {}
                dat['fxx'] = tmp['outp'][0]['fxx'][0]
                dat['pxx'] = tmp['outp'][0]['pxx'][0]
            except:
                print("Error: File not found!")
                continue

            print('Processing S%d B%d M%d ...' % (isubj,iblock,m))


            freqs = np.squeeze(dat['fxx'])
            aper = np.empty([2,dat['pxx'].shape[1]])
            only_gauss = np.zeros([dat['pxx'].shape[0],dat['pxx'].shape[1]])
            full_gauss = np.zeros([dat['pxx'].shape[0],dat['pxx'].shape[1]])


            fm = FOOOFGroup(peak_width_limits=[1, 8], min_peak_height=0.05, max_n_peaks=6)                
            fm._maxfev = 30000
            freq_range = [3, 40]
            fm.fit(freqs, np.transpose(dat['pxx']), freq_range)
            tmp = fm.get_results()

            for isens in range(0,dat['pxx'].shape[1]):
              aper[:,isens] = tmp[isens].aperiodic_params

            F = fm.freqs

            for isens in range(0,dat['pxx'].shape[1]):
              for i in range(0,len(tmp[isens].gaussian_params)):
                  c = tmp[isens].gaussian_params[i][0]
                  w = tmp[isens].gaussian_params[i][2]
                  a = tmp[isens].gaussian_params[i][1]
                  only_gauss[:,isens] = only_gauss[:,isens] + a * np.exp ((-(F-c)**2)/(2*pow(w,2))) 
              b = tmp[isens].aperiodic_params[0]
              e = tmp[isens].aperiodic_params[1]
              full_gauss[:,isens] = only_gauss[:,isens] + b - np.log10(pow(F,e)) 

            scipy.io.savemat('/home/tpfeffer/pupmod/proc/src/pupmod_task_fooof_s%d_m%d_b%d_v%d.mat' % (isubj,m,iblock,v), {'only_gauss': only_gauss, 'full_gauss': full_gauss,  'aper':  aper})

    






        