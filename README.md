# CalPE: compute the Parasitic Energy for CO<sub>2</sub> capture

### References
* [Evaluating different classes of porous materials for carbon capture (2014)](http://xlink.rsc.org/?DOI=C4EE02636E)
* [In silico screening of carbon-capture materials (2012)](http://dx.doi.org/10.1038/nmat3336)

### Authors
* Johanna M. Huck
* Li-Chiang Lin
* Cory M. Simon
* Adam Berger
* Daniele Ongari

### Dependencies

* [pyIAST](https://github.com/CorySimon/pyIAST)
* numpy
* pandas

### Input and run

```
calPE.py Mg-MOF74 914.88 coal -cp 896 -process TPSA -datapath ./test/
```

See the `--help` for the input description.

Use `--log` for printing the debug log file.

#### NB:

* The isotherm data should be put in the `{datapath}/{structure_name}` folder.

* The temperature at which the isotherm data is calculated is automatically
read from the filename `{datapath}/{structure_name}/{adsorbate_name}/{temperature}.csv`.

* Isotherms are fitted using [`pyiast.InterpolatorIsotherm`](https://pyiast.readthedocs.io/en/latest/#interpolatorisotherm)
with `fill_value = uptake.max()`. Therefore, the isotherm should be well
saturated, because for higher pressures the loading is extrapolated as the
maximum uptake.

* The heat of adsorption (HoA) needs to be provided in kJ/mol for all the
loading pressures of the isotherm. It is needed to shift the original isotherm
to a new temperature using the Classius-Clapyeron equation. Note that the HoA
is defined here with a NEGATIVE value.

### Output

In the output, the program prints:

* name of the structure
* parasitic energy (kJ/kg)
* optimal desorption pressure (bar)
* optimal desorption temperature (K)
* fraction of electricity loss (-)
* heat requirement (kJ/kg)
* compression work (kJ/kg)
* mass of CO<sub>2</sub> produced (kg)
* working capacity (mol/kg)
* fraction of CO<sub>2</sub> purity (-)

or, in case of negative working capacity:

* name of the structure
* "Unfeasible process"
