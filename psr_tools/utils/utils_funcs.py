from scipy import integrate
from scipy import interpolate
from scipy import special
from scipy import stats
import numpy as np

def H1(x):
    return -0.97947838884478688 * x - 0.83333239129525072 * x**(1/2) + 0.15541796026816246 * x**(1/3)

def H2(x):
    return -4.69247165562628882e-2 * x - 0.70055018056462881 * x**(1/2) + 1.03876297841949544e-2 * x**(1/3)

def F_curavture(x):
    return 2.149528241534478636710291262 * x**(1/3) * np.exp(H1(x)) + 1.25331413731550025120788 * np.sqrt(x) * np.exp(-x) * (1 - np.exp(H2(x))) 

def F_curavture_exact(x):
    return x * integrate.quad(lambda t: special.kv(5/3, t), x, np.inf)[0]

def is_float(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
  

def inverse_sample_function(dist, Npnts, x_min=-100, x_max=100, n=1e6, random_state=None, **kwargs):
  
  x = np.linspace(x_min, x_max, int(n))
  cumulative = np.cumsum(dist(x, **kwargs)) 
  cumulative -= cumulative.min()
  f = interpolate.interp1d(cumulative/cumulative.max(), x)
  if random_state is not None:
    return f(random_state.random(int(Npnts)))
  else:
    return f(np.random.random(int(Npnts)))

def calculate_width(chi, beta, rho):
    dzeta = beta + chi
    if np.sqrt(np.sin((rho + beta)/2) * np.sin((rho - beta)/2) / np.sin(chi) / np.sin(dzeta)) >= 1:
        return 2 * np.pi
    return 4 * np.arcsin(np.sqrt(np.sin((rho + beta)/2) * np.sin((rho - beta)/2) \
    / np.sin(chi) / np.sin(dzeta)))

def smooth_angle_array(angles, module=180):
    if (module) != 180 and (module != 360):
        print("ARE YOU SHURE THAT module VALUE IS NOT 180 OR 360 degree ???")
    angles_to_process = np.array(angles)
    new_angles = np.zeros_like(angles_to_process)
    new_angles[0] = angles_to_process[0]
    for index in range(1, len(angles_to_process)):
        angle1 = (new_angles[index - 1] // module * module) + angles_to_process[index] % module
        angle2 = (new_angles[index - 1] // module * module) - module + angles_to_process[index] % module
        angle3 = (new_angles[index - 1] // module * module) + module + angles_to_process[index] % module
        if (np.abs(angle1 - new_angles[index - 1]) < np.abs(angle2 - new_angles[index - 1])) and (np.abs(angle1 - new_angles[index - 1]) < np.abs(angle3 - new_angles[index - 1])):
            new_angles[index] = angle1
        elif (np.abs(angle2 - new_angles[index - 1]) < np.abs(angle1 - new_angles[index - 1])) and (np.abs(angle2 - new_angles[index - 1]) < np.abs(angle3 - new_angles[index - 1])):
            new_angles[index] = angle2
        else:
            new_angles[index] = angle3
    return new_angles    

def shift_angle(angle):
    angle = np.asarray(angle)
    scalar_in = (angle.ndim==0)
    angle_vals = np.array(angle, copy=True, ndmin=1, dtype=float)
    for i in range(angle_vals.shape[0]):
        while angle_vals[i] > 90:
            angle_vals[i] -= 180
        while angle_vals[i] < -90:
            angle_vals[i] += 180
    if scalar_in:
        return angle_vals[0]
    else:
        return angle_vals
    
def find_power_law_ind(xs, ys):
    mask = ((np.isnan(xs))+(np.isnan(ys)))
    log_xs = np.log(xs[~mask])
    log_ys = np.log(ys[~mask])
    res = stats.linregress(log_xs, log_ys)
    return res.slope, res.rvalue