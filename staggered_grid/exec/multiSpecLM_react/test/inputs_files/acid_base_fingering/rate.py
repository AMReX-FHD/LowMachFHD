import math

# cgs units
# A = HCl, B = NaOH

# from de Wit's paper
DA = 3.336e-5
DB = 2.129e-5
DAB = DA+DB
print "DAB = %e" % DAB

# from the diameter of water = 3e-8
RAB = 5e-8
print "RAB = %e" % RAB

# Smoluchowski rate (lambda infinite)
kSmol = 4*math.pi*DAB*RAB
print "kSmol = %e" % kSmol

# continuous-space rate (lambda finite)
lam = 1e11
r = lam*RAB**2/DAB
kcont = kSmol*(1-math.tanh(math.sqrt(r))/math.sqrt(r))
print "kcont = %e (r = %e)" % (kcont,r)

# RDME rate
beta = 0.25272
dx = 6.25e-3 
tmp = DAB*dx/beta
kRDME = 1/(1/kcont-1/tmp)
print "kRDME = %e (DAB*dx/beta = %e)" % (kRDME,tmp)

# my rate
kMyrate = 1e-20
print "kMyrate = %e" % kMyrate

# stability limit from reaction
mA = 6.06e-23
mB = 6.64e-23
rhoA0 = 0.0365
rhoB0 = 0.02
nA0 = rhoA0/mA
nB0 = rhoB0/mB
dtReact = 1/kMyrate/max(nA0,nB0)
print "dtReact = %e" % dtReact
