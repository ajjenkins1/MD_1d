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
# harmonic osc
#  v = 0.5*k*(x-xeq)**2
# Morse Potential
  v = (d * (np.exp( -2.0*a*(x-xeq) ) - 2.0*(np.exp(-a*(x-xeq)) )) + c)
  return v

def kinen(v):
  ke = 0.5*mass*v**2
  return ke

def accel(x):
# harmonic osc
#  ac = -k*(x-xeq)/mass
# Morse Potential
  deriv = 2.0*d*a*( np.exp( -2.0*a*(x-xeq) ) )*( (np.exp(a*(x-xeq))) - 1.0)
  ac = -deriv/mass
  return ac 

## Prop methods

def euler(x,v,delT):
  v_new = v + delT*accel(x)
  x_new = x + v*delT
  return (v_new, x_new) 

def velverlet(x,v,delT):
  a = accel(x)
  x_new = x + v*delT + 0.5*a*(delT**2)
  a_new = accel(x_new)
  v_new = v + 0.5*(a + a_new)*delT
  return (v_new, x_new)

######## MAIN PROGRAM ########

if __name__ == '__main__':

# Set some defaults for timestep and trajectory length (overwritten by those specified in input.dat)
  k = 1.e-2
  NStep = 100000 # number of steps to compute
  delT = 1.0 # in a.u.
  dospec = False
  Verlet = False

# First, print a timestamp, set some parameters from the input file
  print('Start-time '+str(strftime("%Y-%m-%d %H:%M:%S")))
  file = open("input_hcl.dat")
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
#   choice of integrator
    if Verlet:
      v,x = velverlet(x,v,delT)
    else:
      v,x = euler(x,v,delT)

    Xlist[istep] = x  #generates a coordinate list#
    Vlist[istep] = v  #generates a velocity list# 
    Elist[istep] = toten(x,v) #generates a total energy list#
  print('Finished trajectory at '+str(strftime("%Y-%m-%d %H:%M:%S"))+'.  Writing results to file')

### OUTPUT ###
# Write results to files for the current trajectory.
  
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
    w = np.real(fftfreq(n,d=delT)*2.0*np.pi*219474.63) # energy list in cm-1
#    w = np.real(fftfreq(n,d=delT)*2.0*np.pi)          # energy list in au 
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


