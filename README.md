
# psr_tools


The projects aims to combine various models and routines for computing pulsar secondary plasma characterisitcs and mean emission profiles into one self-contained Python library.


## Documentation

[Documentation](https://linktodocumentation) - to be done


## Installation

### With pip  

```bash
to be done
```

### From source 
To use pure python moudeles one should just clone this repository and install appropriate python dependencies.  
In oder to build modules with c++ extensions (currently it is only psr_tools.radiative_transfer module) one should install c++ dependencies and run cmake in psr_tools/radiative_transfer folder.
```bash
cd psr_tools
cd radiative_transfer
cmake .
make
```
## Dependencies
Radiative transfer moudle depends on [Eigen3](https://libeigen.gitlab.io/) and [pybind11](https://pybind11.readthedocs.io/en/stable/index.html) c++ libraries.
## Usage/Examples

### Simple profile computation

```python
import numpy as np
import psr_tools
# Create radio pulsar object with observational parameters
PSR = psr_tools.ObservedRadioPulsar('Pulsar', P=0.56, Pdot=1e-15, B12=1.23, chi_deg=26, beta_deg=2.4, freq=1250) 
# Create magnetosphere and emission model
model = psr_tools.radiative_trasfer.FixedHeightModel(multiplicity=3066, gamma=103, Rem=44)
# Create profile calculator
profile_calculator = psr_tools.radiative_trasfer.ProfileCalculator(PSR, model)
# Compute profile
profile = profile_calculator.calculate_profile(-15, 15, 1, 0)
# Plot profile
profile.plot_profile(zoom=True)

```

## Module description
- `plasma_production` - contains various models of secondary plasma production above radio pulsar polar cap
- `data_processing` - contains routins to process observational or modelled data
- `radiative transfer` - contains models of polarisation transfer in radio pulsar magnetophsere. Currently the module is heavily based on [Loki-pulsar-propagation](https://github.com/haykh/Loki-pulsar-propagation) repository



## Feedback

If you have any feedback, please feel free to write to one of the corresponding authors:   
Istomin Arsseniy - istomin15arseniy@gmail.com  
Kniazev Fedor - kniazev.fa@phystech.edu



## License

[MIT](https://choosealicense.com/licenses/mit/)


## Authors

- [@IstominArseniy](https://github.com/IstominArseniy)
- [@kfash](https://github.com/kfash)


