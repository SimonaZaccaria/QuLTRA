import numpy as np
from scipy.optimize import newton, brentq, root_scalar


def zero_algo(admittance, starting_point, final_point):
    step=0.1
    zero_points=[]
    def f(z):
        det_Y,n,_= admittance(2j*np.pi*1e9*z)
        if n % 4 == 0 or n % 4 == 2:
            # n pari: prendi la parte reale
            det = det_Y.real
        else:
            # n dispari: prendi la parte immaginaria
            det= det_Y.imag
        return det
    
    for a in np.arange(starting_point,final_point,step):
        b=a+2*step
        if f(a)*f(b)<0:
            zero=brentq(f,a,b)
            #print(zero)
            det_Y,_,k= admittance(2j*np.pi*1e9*zero)
            if k<=1 and abs(det_Y) < 1e-6:
                zero_points.append(zero)
   
    if not zero_points:
        raise ValueError("No zeros found in the specified interval.")
    
    #R=cxroots.Rectangle([-0.5,0.5],[starting_point, final_point])
    #sol=R.roots(f)
    #zero_points=[zero.imag for zero in sol.roots]
    #zero_points.sort()
    zero_points_no_duplicate=remove_duplicates(zero_points)
    
    return zero_points_no_duplicate

def zero_algo_complete(admittance,guess_points):
    zero_points=[]
    complex_guess_points=[1j*1e9*x0 for x0 in guess_points]

    def f_with_loss(z):
        det_Y,_,_=admittance(2*np.pi*z)
        return det_Y
    
    for x0 in complex_guess_points:
        zero_point=newton(f_with_loss,x0,tol=1.48e-07, maxiter=150)
        if zero_point.real<0: #exclude non stable points
            zero_points.append([zero_point.imag/1e9,zero_point.real/1e6])
    
    return zero_points

def remove_duplicates(values, tol=1e-10):
    values = sorted(values)
    result = []
    for v in values:
        if not result or abs(v - result[-1]) > tol:
            result.append(v)
    return result