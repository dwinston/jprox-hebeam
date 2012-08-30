import re

f = open('sim_COLLISON.txt','r')
fp = open('sim_COLLISION-fancy.txt','w')

# pass annoying non-data line tha starts with \xb3
for line in f:
    if line[0]=='\xb3': break

for line in f:
    if line[0]=='\xb3':
        m = re.split('\xb3',line)
        ionNum, energy, z, y, x, SP = m[1], m[2], m[3], m[4], m[5], m[6]
        fp.write(ionNum+" "+energy+" "+z+" "+y+" "+x+" "+SP+"\n")

f.close()
fp.close()
