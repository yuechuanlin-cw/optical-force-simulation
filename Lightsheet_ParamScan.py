"""
    Parallel computing version of force calculation 
    This script is used for simulation of localized approxiamtion GLMT 
    for calculation of on-axis optical forces for both focused Gaussian beam and light-sheet beam
    Created by Yuechuan Lin 
    Cornell University 
"""
# light-sheet beam: focused Gaussian beam with one axis extended and collimated while the other axis focused. 
# beam radius/waist here indicates the 1/e2 radius of laser intensity beam profile 
import numpy as np 
import matplotlib.pyplot as plt 
import miepython as mp 
import PyMieScatt as pym 
from ILA_BSC import ILA_BSC_func
from Cprz import Cprz_func
import matplotlib.pyplot as plt  
import time  
import xlsxwriter 
import os 
import multiprocessing as mp
import itertools

current_path = os.path.dirname(os.path.abspath(__file__))
StorageFile = current_path+'/Data_LightSheet_ParamScan_ExperimentSimi_17umCOOHMF.xlsx' 

# preparing excel sheet for data storage 
if os.path.exists(StorageFile):
    os.remove(StorageFile)

workbook = xlsxwriter.Workbook(StorageFile)
worksheet_force = workbook.add_worksheet('Force at center of beam')
worksheet_ab = workbook.add_worksheet('AB') 


start_time = time.time()
# define all global parameters 
i = complex(0.,1.)  # define complex unit
P0 = 1 # optical power as 1 Watts 
nsp = 1.68- i*0  # RI of beads 
nenv = 1.34 # RI of medium
mRI = nsp/nenv # dimensionaless RI of sphere relative to ambient medium  
lam = 0.789e-6 # wavelength of laser beam 
k = 2*np.pi/lam # wavenumber in vaccum 
#a = 2.5e-6 # radius of sphere; sometimes, it is determined by the following size parameter x 
#b = np.asarray((1.9e-6,)) # as diameter of sphere
a = np.asarray((1.7e-6,)) # as radius of sphere, choose either a or b 

x = k*a*nenv
#x =  29.753 #k*a*nenv # dimensionaless sphere size parameters 
#print('Particle diameter is a = {0:8.4f} um'.format(x/(k*nenv)*1e6))  # print particles radius
#w0x = 20e-6
#w0y = 20e-6  # beam waist 
#w0 = np.linspace(0.1e-6,8e-6,70) 
w0x = 1.19e-6     # short axis
w0y =  64.55e-6 #58.8e-6 # long axis 

x0 =  0 #np.linspace(-3e-6,3e-6,7)
y0 =  np.linspace(-100e-6,100e-6,480)
z0 =  0 #np.linspace(-10e-6,10e-6,10)

worksheet_force.write(0,0,'ain index')
worksheet_force.write(0,1,'win index')
worksheet_force.write(0,2,'a (um)')
worksheet_force.write(0,3,'w0x (um)')
worksheet_force.write(0,4,'w0y (um)')
worksheet_force.write(0,5,'x0 (um)')
worksheet_force.write(0,6,'y0 (um)')
worksheet_force.write(0,7,'z0 (um)')

worksheet_force.write(0,8,'Cprz (m2)')
worksheet_force.write(0,9,'Force (pN)')

Cprz_v = np.zeros((len(a),len(y0))) 
Force = np.zeros((len(a),len(y0)))
ain_m = np.arange(0,len(a),1)
win_m = np.arange(0,len(y0),1) 


def Cprz_v_fun(xw):
    xw = np.asarray(xw)     
    ain = xw[0]
    win = xw[1]
    print('Now calculating: (a:{0:3d} of {1:3d},w:{2:3d} of {3:3d})'.format(ain+1,len(a),win+1,len(y0)))
    an,bn = pym.Mie_ab(mRI,x[ain])
    w0x_c = w0x #w0[win] 
    w0y_c = w0y #w0[win] 
    x0_c = x0
    y0_c = y0[win]
    z0_c = z0
    Cprz_v[ain,win] = Cprz_func(an,bn,k,[w0x_c,w0y_c],[x0_c,y0_c,z0_c],worksheet=None,NLimit = None)
    Force[ain,win]=(2*P0/(3.0e8*np.pi*(w0x_c*w0y_c)))*(1e12)*Cprz_v[ain,win] # change the unit to be pN
    return [ain,win,a[ain],w0x_c,w0y_c,x0_c,y0_c,z0_c,Cprz_v[ain,win],Force[ain,win]]

if __name__ == "__main__": 
    proces = 20 # mp.cpu_count(); never exhausting the CPU ! maximum is 40 
    chunksize = int(len(a)*len(y0)/proces)
    paramlist = list(itertools.product(ain_m,win_m))
    paramlist_array = np.asarray(paramlist) 
    with mp.Pool(processes=proces) as pool:
        results = pool.map(Cprz_v_fun,paramlist, chunksize)
    
    pool.close()
    pool.join()
    print(np.shape(results))
    print(type(results))

    for list_re in range(np.shape(paramlist_array)[0]):
        tmp_result = results[list_re]
        tmp_result = np.asarray(tmp_result)
        tmp_index = paramlist_array[list_re] 
        ain = tmp_index[0]
        win = tmp_index[1]  
    #    Force[ain,win] = tmp_result[9] 
        worksheet_force.write(ain*len(y0)+win+1,0,tmp_result[0])  # ain index 
        worksheet_force.write(ain*len(y0)+win+1,1,tmp_result[1])   # win index 
        worksheet_force.write(ain*len(y0)+win+1,2,1e6*tmp_result[2]) # a
        worksheet_force.write(ain*len(y0)+win+1,3,1e6*tmp_result[3]) # w0x 
        worksheet_force.write(ain*len(y0)+win+1,4,1e6*tmp_result[4])  # w0y 
        worksheet_force.write(ain*len(y0)+win+1,5,1e6*tmp_result[5])  # x0 
        worksheet_force.write(ain*len(y0)+win+1,6,1e6*tmp_result[6])  # y0 
        worksheet_force.write(ain*len(y0)+win+1,7,1e6*tmp_result[7])  # z0 
        worksheet_force.write(ain*len(y0)+win+1,8,1e6*tmp_result[8])  # Cpr 
        worksheet_force.write(ain*len(y0)+win+1,9,tmp_result[9])  # Force

  
 
    workbook.close()
   # fig,ax = plt.subplots(1,1)
   # ax.imshow(Force,cmap='jet')

    end_time = time.time() 
    print('Elapsed time is : {0:20.10f}'.format(end_time-start_time))
    

   # plt.show() 


