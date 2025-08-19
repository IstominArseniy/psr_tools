"""
list of physical constants 
"""
c = 3*10**10 # speed of light in cm / s
qe = 4.8/10**10 # elementary charge  (CGS units)
me = 9.1e-28 # electron mass in g
re = 2.8179e-13 # classical electron radius in cm
hbar = 1.055e-27 # reduced plank constant in erg * s
Bcr = me**2 * c**3 / qe / hbar # Schwinger magnetic field in G.
Bcr12 = me**2 * c**3 / qe / hbar / 1e12 # Schwinger magnetic field devided by 10^12 G.
e_lambda_bar = hbar / me / c # electron compton wave-length in cm 