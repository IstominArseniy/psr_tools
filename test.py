import numpy as np
from matplotlib import pyplot as plt

from psr_lib.test_package import test_module
from psr_lib import interface

ATNF_db = interface.PSRDataFrame.from_ATNFdb('Data/psrcat.db')
print(ATNF_db.loc['J0002+6216'])
ATNF_db.get_psr_class('J0002+6216').write_to_file()