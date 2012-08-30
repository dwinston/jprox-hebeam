import re

f = open('sim_COLLISION-fancy.txt','r')

maxz = 0;
for line in f:
    m = re.match(r"(.+?) (.+?) (.+?) (.+?) (.+?) (.+?)",line)
    ionNum, energy, z, y, x, SP = m.group(1,2,3,4,5,6)
    if z > maxz:
        maxz = z
print "max depth is",maxz,"Angstroms"
f.close()
