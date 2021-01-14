# slurpy

Python module to solve the 1D, steady, spherical slurry system outlined in Wong et al.
(in prep) (see also Wong et al. 2018).

## Getting Started

### Prerequisites
- [Python](https://www.python.org/)

### Installing
Conda:
```
conda install -c jnywong nondim-slurry
```

Pip:
```
pip install nondim-slurry
```

Git:

Download the latest version of the repository [here](https://github.com/jnywong/nondim-slurry).

## Example scripts

Sample scripts can be found within the module package `$PATH/lib/python3.8/site-packages/slurpy/scripts`.

### Parameter search

1. Open `scripts/parameter_search.py`

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

3. Run `parameter_search.py`

4. Admire the output:

![](https://raw.githubusercontent.com/jnywong/nondim-slurry/master/slurpy/docs/simple_output.png)

### Sensitivity study

1. Open `scripts/sensitivity.py`

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

3. Run `sensitivity.py`

4. Admire the output:

![hello!](https://raw.githubusercontent.com/jnywong/nondim-slurry/master/slurpy/docs/sensitivity_example.png)

## Links
* [PyPI](https://pypi.org/project/nondim-slurry/)
* [Anaconda Cloud](https://anaconda.org/jnywong/nondim-slurry)

## Authors

* [**Jenny Wong**](https://jnywong.github.io/) - *University of Leeds - Institut de Physique du Globe de Paris - Institut des Sciences de la Terre*
* [**Chris Davies**](https://environment.leeds.ac.uk/see/staff/1225/dr-chris-davies) - *University of Leeds*
* [**Chris Jones**](https://eps.leeds.ac.uk/maths/staff/4042/professor-christopher-jones-) - *University of Leeds*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* ERC SEIC
* Del Duca Foundation
* EPSRC Centre for Doctoral Training in Fluid Dynamics

:tada:
