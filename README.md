# How to calculate the Parasitic Energy

### Installation

To use improved_calPE.py, you will first have to install the pyIAST python package.

Instructions on how to download the program can be found here: http://pyiast.readthedocs.io/en/latest/.

Link to the Github page: https://github.com/CorySimon/pyIAST.


2. To run the program, enter the following command into the command line:

```shell
./improved_calPE.py coreDB.yml ABAVIJ coal IAST IAST
```
* `improved_calPE.py` is program to calculate the parasitic energy
* `coreDB.yml` is the database
* `ABAVIJ` is the name of the structure for which you want to calculate the parasitic energy
* `coal` is the composition of the gas (alt. `NG`)
* `IAST` is the mixed-gas model
* `IAST` is the name of the dictionary key in the database

* `improved_calPE.py` is program to calculate the parasitic energy
* `coreDB.yml` is the database
* `ABAVIJ` is the name of the structure for which you want to calculate the parasitic energy
* `coal` is the composition of the gas (alt. `NG`)
* `IAST` is the mixed-gas model
* `IAST` is the name of the dictionary key in the database

#### NB:

* The isotherm data should be put in the `ccsdata/{structure_name}` folder

* The temperature at which the isotherm data is calculated is automatically read from the
filename `ccsdata/{structure_name}/{adsorbate_name}/{temperature}.csv`

* Make sure to specify the `{structure_name}` in the coreDB.yml file:
```
Structures:
  {structure_name}
    Adsorbates:
```
### Output

The output you get should look like something like this:

```shell
ABAVIJ [914.4198, 2026.5, 333.0, 0.248214, 192.9399, 721.48,43.1682, 0.98971, 0.924254]
```

* `ABAVIJ` is the name of the structure
* `914.4198` is the parasitic energy in kJ/kg
* `2026.5` is the pressure (Pa)
* `333.0` is the final temperature (K)
* `0.248214` is the fraction of final electricity loss (-)
* `192.9399` is the final heat requirement (kJ/kg)
* `721.48` is the final compression work (kJ/kg)
* `43.1682` is the final mass of CO2 produced (kg)
* `0.98971` is the final working capacity (mol/kg)
* `0.924254` is the fraction of the final purity CO2 (-)
