import math
import sys

if (len(sys.argv)!=3):
  print "ERROR: rate_multiplier and D_Fick expected"
  sys.exit()
else:
  RM = float(sys.argv[1]) 
  DFICK = float(sys.argv[2])

hydro_grid_output = "Schlogl.S_k.pair=1.Re.dat"
output = "res.Sk_xi"

# values from Schlogl monomodal
r0 = 1.75
G0 = 3.55
nb = 1.
C = G0/nb/r0

# conversion factor from k to xi
factor = math.sqrt(DFICK/(RM*r0))

out = open(output,"w")

with open(hydro_grid_output) as inf:
  for line in inf:

    linesplit = line.split()

    # if nothing
    if (len(linesplit)==0):	
      continue

    # if starting with '#'
    if (linesplit[0][0]=='#'):
      continue

    val1 = float(linesplit[0])
    val2 = float(linesplit[1])

    xi = factor*val1
    Sk_exact = (C+xi*xi)/(1+xi*xi)

    out.write("%g\t%g\t%g\n"%(xi,val2,Sk_exact))

out.close()
