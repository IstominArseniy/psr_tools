from scipy import integrate
from scipy import special
import numpy as np

def H1(x):
    return -0.97947838884478688 * x - 0.83333239129525072 * x**(1/2) + 0.15541796026816246 * x**(1/3)

def H2(x):
    return -4.69247165562628882e-2 * x - 0.70055018056462881 * x**(1/2) + 1.03876297841949544e-2 * x**(1/3)

def F_curavture(x):
    return 2.149528241534478636710291262 * x**(1/3) * np.exp(H1(x)) + 1.25331413731550025120788 * np.sqrt(x) * np.exp(-x) * (1 - np.exp(H2(x))) 

def F_curavture_exact(x):
    return x * integrate.quad(lambda t: special.kv(5/3, t), x, np.inf)[0]

