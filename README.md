# slurpy

Solve the 1D, steady, spherical slurry system outlined in Wong et al.
(in prep) (see also Wong et al. 2018).

## Getting Started

### Prerequisites
- [Python](https://www.python.org/)

### Installing
Conda:
```
conda install nondim-slurry
```

Pip:
```
pip install nondim-slurry
```

Git:
Download the latest version of the repository [here](https://github.com/jnywong/nondim-slurry).

### A simple example

Sample scripts can be found within the module package `slurpy/scripts`.

1. Open `parameter_search.py`

2. Enter some input parameters. For example, try:

```
# %% MODEL INPUTS
# Show plots?
plotOn=1 # show temp, xi, solid flux and density profiles

# Input parameters
layer_thicknesses=np.array([150e3]) # (m)
thermal_conductivities=np.array([100.]) # (W m^-1 K^-1)
icb_heatfluxes=np.array([3.4]) # (TW)
csb_heatfluxes=np.array([7.4]) # (TW)

h=0.05 # stepsize of heat flux through parameter space
```

![](https://raw.githubusercontent.com/jnywong/nondim-slurry/master/slurpy/docs/simple_example.png)

3. Run `parameter_search.py`

4. Admire the output:

![](https://raw.githubusercontent.com/jnywong/nondim-slurry/master/slurpy/docs/simple_output.png)

### Example: Sensitivity study

1. Open `sensitivity.py`

2. Enter some input parameters. For example, try:

```
# %% MODEL INPUTS
# Save plot?
saveOn=0

# Input parameters
layer_thickness=150e3 # (m)
thermal_conductivity=100. # (W m^-1 K^-1)
icb_heatflux=2.5 # (TW)
csb_heatflux=5.0 # (TW)
h=0.05 # stepsize of heat flux through parameter space

# Sensitivity study
csb_temp = np.arange(4500.,6100.,100) # (K)
csb_oxy = np.arange(2,12.5,0.5) # (mol.%)
sed_con= np.array([1e-5,1e-4,1e-3,1e-2,1e-1]) # (kg s/m^3) pre-factor in sedimentation coefficient, b(phi)
```

![](https://raw.githubusercontent.com/jnywong/nondim-slurry/master/slurpy/docs/sensitivity.png)

3. Run `sensitivity.py`

4. Admire the output:

![hello!](https://raw.githubusercontent.com/jnywong/nondim-slurry/master/slurpy/docs/sensitivity_example.png)

## Links

* [PyPI](https://test.pypi.org/project/slurpy/)

## Authors

* **Jenny Wong** - *Institut de Physique du Globe de Paris*
* **Chris Davies** - *University of Leeds*
* **Chris Jones** - *University of Leeds*

## License

This project is licensed under the MIT License - see the [license.md](LICENSE.md) file for details

## Acknowledgments

* Del Duca Foundation
* EPSRC Centre for Doctoral Training in Fluid Dynamics

:tada:
