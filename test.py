import numpy as np
from matplotlib import pyplot as plt
import tqdm
from psr_lib.test_package import test_module
from psr_lib import interface
import time

psr = interface.RadioPulsar.from_file('test_pulsar.toml')
t1 = time.time()
M = psr.curvature_multipliticy(0.5, 3e7)
t2 = time.time()
print(M, t2-t1)