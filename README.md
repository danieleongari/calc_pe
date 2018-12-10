# Compute the Parasitic Energy for CO2 separation

### Installation

To use improved_calPE.py, you will first have to install the pyIAST python package.

Instructions on how to download the program can be found here: http://pyiast.readthedocs.io/en/latest/.

Link to the Github page: https://github.com/CorySimon/pyIAST.

# Input and run

See the `--help` for the input description.

Use `--log` for printing the debug log file.

#### NB:

* The isotherm data should be put in the `ccsdata/{structure_name}` folder

* The temperature at which the isotherm data is calculated is automatically
read from the filename `ccsdata/{structure_name}/{adsorbate_name}/{temperature}.csv`

### Output

In the output, the program prints:

* the name of the structure
* the parasitic energy (kJ/kg)
* the optimal desorption pressure (bar) #changed from (Pa)!
* the optimal desorption temperature (K)
* the fraction of electricity loss (-)
* the heat requirement (kJ/kg)
* the compression work (kJ/kg)
* the mass of CO2 produced (kg)
* the working capacity (mol/kg)
* the fraction of CO2 purity (-)

#### NB:

* "-1" values stand for "None"
