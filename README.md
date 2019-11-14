# calc_pe

calc_pe computes the parasitic energy for CO<sub>2</sub> capture from CO<sub>2</sub> and N<sub>2</sub>  isotherms.

### Install
```
git clone https://github.com/danieleongari/calc_pe.git
cd calc_pe
pip install .
```

### Input and run

```
$ calc_pe Mg-MOF74 coal -rho 914.88 -cp 896 -process TPSA -datapath ./tests/
```

See `calc_pe --help` for the input description.

Use `calc_pe --log` for printing the debug log file.

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

* You can provide density and cp as single value files `cp.csv` and `rho.csv`:
see the tests as example.

* For testing the minimal inputs are:
```
$ cd tests/
$ calc_pe Mg-MOF74 coal
$ calc_pe HKUST-1 coal
```

### Output

In the output, the program prints:

```
Mg-MOF74: PE(MJ/kg)= 0.867: Pd(bar)= 0.01 Td(K)= 333.0 EL(-) = 0.235 Q(MJ/kg)= 0.124 Wcomp(MJ/kg)= 0.743 WCv(kg/m3)= 114.655 WCg(kg/kg)= 0.193 pur(-)= 0.967
```

* Name of the adsorbent
* `PE(MJ/kg)`: parasitic energy per kg of CO<sub>2</sub> (Note: PE=Q+Wcomp)
* `Pd(bar)`: optimal desorption pressure
* `Td(K)`: optimal desorption temperature
* `EL(J/J)`: fraction of electricity loss
* `Q(MJ/kg)`: heat requirement
* `Wcomp(MJ/kg)`: compression work
* `WCv(kg/m3)`: volumetric working capacity, i.e.,
mass of CO<sub>2</sub> produced per m<sup>2</sup> of bed,
considering `-vf` void fraction.
* `WCg(kg/kg)`: gravimetric working capacity, i.e.,
mass of CO<sub>2</sub> produced per kg of bed,
considering `-vf` void fraction.
* `pur(mol/mol)`: molar fraction of CO<sub>2</sub> final purity (-)

A warning is printed in case of negative working capacity
for all the tested desorption conditions, e.g.:

```
$ calc_pe HKUST-1 air
HKUST-1: Unfeasible process!
```

#### NB:

* The Henry coefficient for CO<sub>2</sub> is a good pre-screening parameter

* The working capacity is also very important, since it allows for less cycles
using the same amount of adsorbent (or less adsorbent needed with the same
cycles).

* The final CO<sub>2</sub> purity is less than the imposed purity, `-yd`
(default: 0.99): we use the `yd` value as an approximation of the gas phase at
desorption to get the uptake in the adsorbent at the desorption condition
(using IAST). Note that the PE is not very sensitive to `yd`
(see [Joos et al. (2016)](http://doi.org/10.1039/c6fd00031b))
and there is not a motivated need for reiteration.
The final CO<sub>2</sub> purity is computed as the working capacity of
CO<sub>2</sub> over the sum of the working capacities of both CO<sub>2</sub>
and N<sub>2</sub>.

* By default the program prints the results for optimal PE (i.e., the lowest).
However, one can search for other optimal parameters by using the `-opt` command:
lowest `Q` if he is not interest in compressing the CO<sub>2</sub>,
highest working capacity (`WC`) or highest CO<sub>2</sub> final purity (`pur`).
*Note that these may not be anymore optimization problems, returning just
the max/min T and P conditions.*


### Dependencies

calc_pe uses:

* [pyIAST](https://github.com/CorySimon/pyIAST)
* numpy
* pandas

### References

If you use calc_pe, please consider citing:

* [Evaluating different classes of porous materials for carbon capture (2014)](http://doi.org/10.1039/C4EE02636E)
* [In silico screening of carbon-capture materials (2012)](http://dx.doi.org/10.1038/nmat3336)

This program has been used in:

* [Building a Consistent and Reproducible Database for Adsorption Evaluation in Covalent–Organic Frameworks](https://pubs.acs.org/doi/abs/10.1021/acscentsci.9b00619)

### Authors
* Johanna M. Huck
* Li-Chiang Lin
* Cory M. Simon
* Adam Berger
* Daniele Ongari (restyling, December 2018)
