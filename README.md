
# psr_tools


The projects aims to combine various models and routines for computing pulsar secondary plasma characterisitcs and mean emission profiles into one self-contained Python library.


## Documentation

[Documentation](https://linktodocumentation) - to be done


## Installation

With pip:  

```bash
to be done
```

From source:  
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

```python
import psr_lib

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


