import numpy as np
from scipy.optimize import root
from scipy.constants import epsilon_0, mu_0

c=299792458 #m/s
e = 1.60217657e-19  # electron charge
h = 6.62606957e-34  # Plank's
hbar=h/2/np.pi
phi0=h/(2*e)
epsilon_r=11.9
epsilon_eff=(1+epsilon_r)/2 #for silicon substrate
v=c/np.sqrt(epsilon_eff)
step=0.1 #step for frequency sweep
newton_tol=1.48e-7 #tolerance for newton method
maxiter=150
k_max=100 #MHz