###### MODULES TO BE USED IN THIS SCRIPT #######
from __future__ import division
import math
import numpy as np
from time import strftime
from scipy.fftpack import fft, fftfreq
import matplotlib.pyplot as plt

###### FUNCTION DEFINITIONS #######

def toten(x,v):
  E = kinen(v) + poten(x)
  return E

def poten(x):
  v = 0.5*k*(x-xeq)**2
  return v

def kinen(v):
  ke = 0.5*mass*v**2
  return ke

def accel(x):
  a=-k*(x-xeq)/mass
  return a 

## Prop methods

def euler(x,v,delT):
  v_new=v+delT*accel(x)
  x_new=x+v*delT
  return (v_new, x_new) 

######## MAIN PROGRAM ########

if __name__ == '__main__':

#Set some defaults for timestep and trajectory length (overridden by those specified in input.dat)
  k = 1.e-2
  NStep=100000 # number of steps to compute
  delT =1.0 # in a.u.
  dospec = False

#First, print a timestamp, set some parameters from the input file
  print('Start-time '+str(strftime("%Y-%m-%d %H:%M:%S")))
  file = open("input.dat")
  for line in file:
    exec(line)

# Initialize some arrays to store the position, velocity, total energy, and time at each step.
  Elist =[0]*Nstep 
  Xlist =[0]*Nstep 
  Vlist =[0]*Nstep
  tlist =[0]*Nstep

# Perform numerical integration of the classical EOM

  x=xinit
  v=vinit
  for istep in range(0,Nstep):
    tlist[istep]=delT*istep
    v,x=euler(x,v,delT)
    Xlist[istep] = x  #generates a coordinate list#
    Vlist[istep] = v  #generates a velocity list# 
    Elist[istep] = toten(x,v) #generates a total energy list#
  print('Finished trajectory at '+str(strftime("%Y-%m-%d %H:%M:%S"))+'.  Writing results to file')

### OUTPUT ###
# Write results to files for the current trajectory.

#  MD_coords = open('MD_coords.xyz', 'wt')
# Currenly assuming H2 as system.  Change accordingly for accurate visualization.
#  for item in Xlist:
#    MD_coords.write("2\n \n H 0.0 0.0 0.0 \n H 0.0 0.0 %.8f\n" % item)
#  MD_coords.close() 
  
  np.savetxt('MD_velocity.txt', np.transpose([tlist,Vlist]))  
  np.savetxt('MD_toten.txt', np.transpose([tlist,Elist]))  
  np.savetxt('MD_x.txt', np.transpose([tlist,Xlist]))  

# Energy Conservation
  Energy_Error = Elist[-1] - Elist[0]  
  print('Energy Conservation Error ', Energy_Error)

  if (dospec):
# FFT the velocity and write spectrum (frequency in first column, amplitude in second) to file
    fw = np.real(abs(fft(Vlist)))
    n = len(Vlist)
#    w = np.real(fftfreq(n,d=delT)*2.0*np.pi*219474.63) # energy list in cm-1
    w = np.real(fftfreq(n,d=delT)*2.0*np.pi) # energy list in au 
    np.savetxt('MD_spectrum.txt', np.transpose([w,fw]))  
  print('Done!!!! '+str(strftime("%Y-%m-%d %H:%M:%S")))
 
### Display plots ###
# Spectrum
  plt.plot(w,fw)
  plt.xlabel("Energy")
  plt.ylabel("Amplitude")
  plt.show()
# x position in time
#  plt.plot(tlist,Xlist)
#  plt.xlabel("Bond Length (Bohr)")
#  plt.ylabel("Step")
#  plt.show()


