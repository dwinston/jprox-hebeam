''' Put the material you want at the end '''

# PMMA -> (C5O2H8)n
names = 'H','C','O'
numbers = 1, 6, 8
weights = 1.0079, 12.0107, 15.999
stoichs = 8, 5, 2
d = 1.19 # total density, in g/cm^3

#HSQ
names = 'H', 'Si', 'O'
numbers = 1, 14, 8
weights = 1.008, 28.09, 16.000 # atomic weights of H, Si, O
stoichs = 8, 8, 12 # stoichiometry: H8Si8O12
d = 1.4 # total density, in g/cm^3

if len(weights) != len(stoichs) and len(weights) != len(names):
    print "error: weights and stoichs must be equal in size"
    import sys; sys.exit(1)
material_weight = sum([s*w for s,w in zip(stoichs,weights)])

# A human-friendly format
for n,w,s in zip(names, weights, stoichs):
    print "density of species",n,"is", d*w*s/ material_weight

# For copy-and-paste into a simulation input file
print len(names)
for Z,w,s in zip(numbers, weights, stoichs):
    print Z
    print w
    print d*w*s / material_weight


