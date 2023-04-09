from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt


def tobohr(x):
    x = x*1.8897161646320724
    return x

def morse(x, d, a, u , c):
    return (d * (np.exp(-2*a*(x-u))-2*np.exp(-a*(x-u))) + c)

def fit_morse():
    # Initial parameters for fit guess. Take some from pes.dat.
    # Minimum energy index
    Minloc = np.argmin(File_data[:,1])
    # Use the min energy and bond length in the guess 
    tstart = [File_data[Minloc,1], 1, File_data[Minloc,0], 0]

    popt, pcov = curve_fit(morse, File_data[:,0], File_data[:,1], p0 = tstart,  maxfev=40000000)
    return popt

# Text file data converted to integer data type
File_data = np.loadtxt("pes.dat", dtype=float)
# Convert to Bohr?
File_data[:,0] = tobohr(File_data[:,0])
popt = fit_morse()

### Print Parameters
# print('Dissociation energy, d:',popt[0])
# print('Curvature parameter, a:',popt[1])
# print('Equilibrium bond length, u:',popt[2])
# print('Shift, c:',popt[3])
# to file:
np.savetxt('morse_params.txt', np.transpose([popt]))

### Test it worked - Plot computed PES together with Morse fitted surface
# Generate Points for Morse
t = np.linspace(1.5,10)
yfit = morse(t,popt[0], popt[1], popt[2], popt[3])

plt.plot(File_data[:,0], File_data[:,1],"ro")
plt.plot(t, yfit)

plt.xlabel("Bond Length (Bohr)")
plt.ylabel("Energy (Hartree)")

plt.show()

### Test end





