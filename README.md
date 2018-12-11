# Compute the Parasitic Energy for CO<sub>2</sub> separation

### References
* [Evaluating different classes of porous materials for carbon capture (2014)](http://xlink.rsc.org/?DOI=C4EE02636E)
* [In silico screening of carbon-capture materials (2012)](http://dx.doi.org/10.1038/nmat3336)

### Authors
* Johanna M. Huck
* Li-Chiang Lin
* Cory M. Simon
* Adam Berger

### Dependencies

* [pyIAST](https://github.com/CorySimon/pyIAST)
* numpy
* pandas

### Input and run

See the `--help` for the input description.

Use `--log` for printing the debug log file.

#### NB:

* The isotherm data should be put in the `ccsdata/{structure_name}` folder.

* The temperature at which the isotherm data is calculated is automatically
read from the filename `ccsdata/{structure_name}/{adsorbate_name}/{temperature}.csv`.

* Isotherms are fitted using [`pyiast.InterpolatorIsotherm`](https://pyiast.readthedocs.io/en/latest/#interpolatorisotherm)
with the max loading as `fill_value`. Therefore, the isotherm should be well
saturated, because for higher pressures the loading is extrapolated as the
maximum loading.

* The heat of adsorption (HoA) needs to be provided in kJ/mol for all the
loading pressures of the isotherm. It is needed to shift the original isotherm
to a new temperature using the Classius-Clapyeron equation. Note that the HoA
is defined here with a positive value.

### Output

In the output, the program prints:

* the name of the structure
* the parasitic energy (kJ/kg)
* the optimal desorption pressure (bar) #changed from (Pa)!
* the optimal desorption temperature (K)
* the fraction of electricity loss (-)
* the heat requirement (kJ/kg)
* the compression work (kJ/kg)
* the mass of CO<sub>2</sub> produced (kg)
* the working capacity (mol/kg)
* the fraction of CO<sub>2</sub> purity (-)

#### NB:

* "-1" values stand for "None"
