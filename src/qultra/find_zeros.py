import numpy as np
from scipy.optimize import newton, brentq
from mpmath import findroot,mpc
import warnings


try:
    from . import constants   # Importa il modulo intero
except ImportError:
    import constants  



def zero_algo(admittance, starting_point, final_point):
    
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
    
    for a in np.arange(starting_point,final_point,constants.step):
        b=a+2*constants.step
        if f(a)*f(b)<0:
            zero=brentq(f,a,b)
            #print(zero)
            det_Y,_,k= admittance(2j*np.pi*1e9*zero)
            if k<=1 and abs(det_Y) < 1e-3:
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
    #zero_points1=[]
    zero_points_final=[]
    complex_guess_points=[1j*1e9*x0 for x0 in guess_points]
    


    def f_with_loss(z):
        det_Y,_,_=admittance(2*np.pi*z)
        return det_Y
    
    
    for x0 in complex_guess_points:
        try:
            zero_point=newton(f_with_loss,x0,tol=constants.newton_tol, maxiter=constants.maxiter)
            
        except RuntimeError as e_newton:
            warnings.warn(f"Newton's method failed for initial guess {x0}: {e_newton}")
            zero=x0
            for i in range(101):
                zero=alternative_zero_finding(admittance,zero,i*0.01)
            zero_point=complex(zero)
            
        if zero_point.real<0: #exclude non stable points
            if -2*zero_point.real<=constants.k_max*1e6: #exclude points with k>k_max
                if zero_point.imag >= 0: #positive frequency
                    zero_points.append([zero_point.imag/1e9,zero_point.real/1e6])
            else:
                warnings.warn(f"Zero point {zero_point} has k > k_max, excluding it.")   


    """
    if there_are_duplicates(zero_points, tol=1e-10):
        
        for x0 in complex_guess_points:
            try:
                zero_point1 =  newton(f_with_loss,x0-2*np.pi*50e6,tol=constants.newton_tol, maxiter=constants.maxiter)

            except RuntimeError:
                zero1=x0-2*np.pi*50e6
                for i in range(100):
                    zero1=newton2(admittance,zero1,(i+1)*0.01)
                zero_point1=zero1

            if zero_point1.real<0: #exclude non stable points
                zero_points1.append([zero_point1.imag/1e9,zero_point1.real/1e6])
    """           
    zero_points_final=zero_points#+zero_points1
    zero_points_final_reduced=remove_complex_duplicates(zero_points_final)
    
    
    return zero_points_final_reduced

def remove_duplicates(values, tol=1e-10):
    values = sorted(values)
    result = []
    for v in values:
        if not result or abs(v - result[-1]) > tol:
            result.append(v)
    return result

def remove_complex_duplicates(values, tol=1e-10):
    values = sorted(values, key=lambda x: x[0])  # Sort by the real part
    result = []
    for v in values:
        if not result or abs(v[0] - result[-1][0]) > tol:
            result.append(v)
    return result

def there_are_duplicates(values, tol=1e-10):
    values = sorted(values, key=lambda x: x[0])
    for i in range(1, len(values)):
        if abs(values[i][0] - values[i - 1][0]) <= tol:
            return True
    return False



def alternative_zero_finding(admittance, guess_point, weight):
    def f(z):
        # Cast immediato così NumPy vede solo complex Python
        z_c = complex(z)

        det_Y, n, _ = admittance(2*np.pi*z_c)
        
        if n % 4 == 0 or n % 4 == 2:
            val = det_Y.real + 1j * weight * det_Y.imag
        else:
            val = 1j * det_Y.imag + weight * det_Y.real
        
        # Restituisco mpc per mpmath
        return mpc(val)

    zero_point = findroot(f, mpc(complex(guess_point)), solver='muller')
    return zero_point
