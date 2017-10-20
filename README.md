## Polfit
This will fit a polynomial to the potential energy of a diatomic as a function of separation,
then perform a simple Dunham analysis on it. 

### Dependencies

- Python 3.5 or higher
- Numpy
- texttable
- matplotlib

## Running the script
The input is simply a tab-delimited text file of data, which must include a header row with column names.
See for example `cl2test.dat`. 

The script is then run as (assuming coordinates are in BOHR):

`python3 polfit.py -f [input file] -mu [reduced mass of diatomic in atomic mass units]`

If your coordinates are in ANGSTROM instead, pass the flag

`-angstrom True`

If you want to plot one or more of the polynomial fits, pass

`-plot [list of col. numbers]`

Finally, you can change the order of polynomial fitted with

`-order [order of poly]`

Note that this must be >= 6 for the anharmonicity constants to be calculated.

e.g. for the `cl2test.dat` example (which is in Angstrom):

`python3 polfit.py -f cl2test.dat -mu 17.7265 -angstrom True`

would give the usual analysis as per MOLPRO, while:

`python3 polfit.py -f cl2test.dat -mu 17.7265 -angstrom True -plot 1 2 -order 5`

would do the analysis with a 5th-order polynomial, and plot both the EHF and ECC
(columns 1 and 2 of the input). Changing the plot flag to `-plot 2` would just plot
the ECC fit. 

