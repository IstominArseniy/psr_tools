import numpy as np
from matplotlib import pyplot as plt

from psr_lib.test_package import test_module
from psr_lib import interface

ATNF= interface.ATNFDatabase.from_ATNFdb('Data/psrcat.db')
print(ATNF.loc['J0002+6216'])