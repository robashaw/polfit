import argparse
import numpy as np
import matplotlib.pyplot as plt
from texttable import Texttable

# Conversion factors
TO_CM = 219474.63067
TO_EV = 27.2113839
TO_BOHR = 1.88973
TO_ANGSTROM = 0.5291761
FORCE_MASS = 1822.88853

def polfit(x, y, n=6):
    "Fits a polynomial of order n to the set of (x [Bohr], y [Hartree]) coordinates given, and calculates data necessary for a Dunham analysis" 
    # Find best guess at minimum and shift coordinates
    xref = x[np.argmin(y)]
    xshift = x - xref
    
    # Fit polynomial to shifted system
    z = np.polyfit(xshift, y, n)
    p = np.poly1d(z)
    
    # Find the true minimum by interpolation, if possible
    xmin = min(xshift)-0.1
    xmax = max(xshift)+0.1
    crit_points = [x.real for x in p.deriv().r if np.abs(x.imag) < 1e-8 and xmin < x.real < xmax] 
    if len(crit_points) == 0:
        print("MINIMUM NOT FOUND")
        # Set outputs to default values
        re = xref
        pt = [0.0]*7
        
    else:
        dx = crit_points[0]
        re = xref + dx # Equilibrium geometry
        
        # Calculate 0th - 6th Taylor series coefficients at true minimum
        pt = []
        for i in range(7):
            pi = p.deriv(i)(dx)/np.math.factorial(i)
            pt.append(pi)
    
    # Return fitted polynomial, x-shift, equilibrium bond length,
    # and Taylor series coefficients
    return p, xref, re, pt
            
def dunham(re, pt, mu, n=6, Emax = 0):
    "Performs a Dunham analysis on a diatomic, given equilibrium geometry Re, the first 7 Taylor series coefficients pt, and the reduced mass mu"
    An = mu * FORCE_MASS
    
    # Energy at minimum, first rotational constant, and first vibrational constant
    Ee = pt[0]
    Be = 0.5 * TO_CM / (An * re**2)
    We = TO_CM * np.sqrt(2.0 * np.abs(pt[2]) / An)
    
    # Compute normalised derivatives
    npt = []
    for i in range(4):
        npi = (pt[i+3]/pt[2])*re**(i+1)
        npt.append(npi)
    
    # Second rotational constant
    Ae = -6.0 * Be**2 * (1.0 + npt[0]) / We
    
    # First anharmonic corrections require n >= 6
    Wexe = 0.0
    Weye = 0.0
    if n > 5: 
        Wexe = -1.5 * (npt[1] - 1.25*npt[0]**2) * Be
        Weye = 0.5 * (10.0*npt[3] - 35.0*npt[0]*npt[2] - 8.5*npt[1]**2 + 56.125*npt[1]*npt[0]**2 - 22.03125*npt[0]**4)*Be**2/We
    
    # Dissociation energies
    De = 0.0
    D0 = 0.0 
    if Emax != 0:
        De = (Emax - Ee) * TO_EV
        D0 = De - 0.5 * (We - 0.5*Wexe) * TO_EV/TO_CM
    
    return Ee, Be, Ae, We, Wexe, Weye, De, D0
    
def plot_polfit(x, ys, ps, xshifts, names):
    "Plots the results of one or more polfits"
    xp = np.linspace(min(x)-0.05, max(x)+0.05, 100)

    # Plot all curves
    for i in range(len(ys)):
         plt.plot(x, ys[i], '.')
         plt.plot(xp, ps[i](xp-xshifts[i]), '-', label=names[i])
    
    plt.xlabel('R (Bohr)')
    plt.ylabel('Energy (Ha)')
    plt.legend()
    plt.show()

def read_input_file(filename):
    "Reads a multi-column file of data"
    names = []
    x = []
    ys = []

    with open(filename, 'r') as f: 
        first_line = f.readline()
        names = [str(v) for v in first_line.split()]
        ncols = len(names)
        if ncols > 1:
            for y in range(ncols-1):
                ys.append([])
            for line in f:
                cols = [float(v) for v in line.split()]
                if len(cols) >= ncols: 
                    x.append(cols[0])
                    for y in range(ncols-1):
                        ys[y].append(cols[y+1])
        else:
            print("Input file does not have enough columns!")
        
    return names, x, ys


parser = argparse.ArgumentParser()
parser.add_argument("-f", dest="filename", required=True, type=str, help="File with table of data")
parser.add_argument("-mu", dest="mu", required=True, type=float, help="The reduced mass of the diatomic, in atomic mass units")
parser.add_argument("-plot", dest="plot_ix", required=False, type=int, nargs="+", default=[], help="Plot the resulting polynomials for these columns")
parser.add_argument("-order", dest="poly_order", required=False, type=int, default=6, help="Order of the polynomial to fit (>=6 required for Dexe and Deye)")
parser.add_argument("-angstrom", dest="angstrom", required=False, type=bool, default=False, help="True if distances are in Angstrom (default is Bohr)")
args = parser.parse_args()

names, x, ys = read_input_file(args.filename) 
x = np.array(x)
for y in ys:
    y = np.array(y)

if args.angstrom is True:
    x *= TO_BOHR

ps = []     # Polynomial fits
xrefs = []  # Reference minima for data
results = []

i = 1
for y in ys:
    p, xref, re, pt = polfit(x, y, args.poly_order)
    ps.append(p)
    xrefs.append(xref)
    
    E, B, A, W, Wx, Wy, De, D0 = dunham(re, pt, args.mu)
    results.append([names[i], E, re, re*TO_ANGSTROM, B, A, W, Wx, Wy, De, D0])
    i += 1

table = Texttable()
headings = ['Name', 'Emin (Ha)', 'Re (Bohr)', 'Re (Ang)', 'B (cm-1)', 'A (cm-1)', 'We (cm-1)', 'Wexe', 'Weye', 'De', 'D0']
table.header(headings)
table.set_cols_dtype(['t', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f'])
table.set_precision(6)
table.set_cols_width([5, 11, 8, 8, 8, 8, 8, 8, 8, 8, 8])

for r in results:
    table.add_row(r)
print(table.draw())


if len(args.plot_ix) > 0:
    ysubset = []
    psubset = []
    namesubset = []
    xshifts = []
    for ix in args.plot_ix:
        if (ix-1) < len(ys):
            ysubset.append(ys[ix-1])
            psubset.append(ps[ix-1])
            namesubset.append(names[ix])
            xshifts.append(xrefs[ix-1])
    
    plot_polfit(x, ysubset, psubset, xshifts, namesubset)


        
    

