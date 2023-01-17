#一次处理一个文档
#包含所有的configuration

#1_rearrange
#load file
fread1 = open('long.lammpstrj', 'r')
longread1 = fread1.readlines()
fread2 = open('short.lammpstrj', 'r')
shortread1 = fread2.readlines()
# strip and split data
longread2 = []
shortread2 = []
for i in longread1:
    first = i.strip('/n')
    second = first.split()
    longread2.append(second)
for i in shortread1:
    first = i.strip('/n')
    second = first.split()
    shortread2.append(second)

#delete title, 8 means the length of one line in dump file
longdump = []
shortdump = []
for i in range(0, len(longread2), 1):
    if len(longread2[i]) == 8:
        longdump.append(longread2[i])
for i in range(0, len(shortread2), 1):
    if len(shortread2[i]) == 8:
        shortdump.append(shortread2[i])

#add chain number, the long and short chain length should be changed
longchainnumber = []
shortchainnumber = []
longchainlength = 120
shortchainlength = 40
totalparticleslong = 6240 #atoms number
totalparticlesshort = 6240  #atoms number
for i in range(1, len(longdump)//totalparticleslong+1, 1):
    p = 0
    for j in range(1, totalparticleslong//longchainlength+1, 1):
        p = p + 1
        for k in range(1, longchainlength+1, 1):
            longchainnumber.append(p)
for i in range(1, len(shortdump)//totalparticlesshort+1, 1):
    p = 0
    for j in range(1, totalparticlesshort//shortchainlength+1, 1):
        p = p + 1
        for k in range(1, shortchainlength+1, 1):
            shortchainnumber.append(p)

#write into a new file
fwrite1 = open('long_rearrange.txt', 'w')
fwrite2 = open('short_rearrange.txt', 'w')
a = 'chainN ID TYPE x y z ix iy iz ' + '\n'
d = 'chainN ID TYPE x y z ix iy iz ' + '\n'
fwrite1.write(a)
fwrite2.write(d)
for i in range(0, len(longchainnumber), 1):
    b = str(longchainnumber[i]).replace('[', '').replace(']', '')
    b = b.replace("'", '').replace(',', '') + ' '
    c = str(longdump[i]).replace('[','').replace(']','')
    c = c.replace("'",'').replace(',','') + '\n'
    fwrite1.write(b)
    fwrite1.write(c)
for i in range(0, len(shortchainnumber), 1):
    e = str(shortchainnumber[i]).replace('[', '').replace(']', '')
    e = e.replace("'", '').replace(',', '') + ' '
    f = str(shortdump[i]).replace('[','').replace(']','')
    f = f.replace("'",'').replace(',','') + '\n'
    fwrite2.write(e)
    fwrite2.write(f)
fwrite1.close()
fwrite2.close()

#2_pandas_particles_in_pore
import pandas as pd

#load files
dfl = pd.read_table('long_rearrange.txt', delim_whitespace=True)
dfs = pd.read_table('short_rearrange.txt', delim_whitespace=True)

#calculate real coordination
xbox = 50
ybox = 50
zbox = 101
dfl['realx'] = dfl.x + xbox*dfl.ix
dfl['realy'] = dfl.y + ybox*dfl.iy
dfl['realz'] = dfl.z + zbox*dfl.iz
dfs['realx'] = dfs.x + xbox*dfs.ix
dfs['realy'] = dfs.y + ybox*dfs.iy
dfs['realz'] = dfs.z + zbox*dfs.iz

#select particles in pore
dfl_in_pore = dfl.loc[dfl['realz'] < 0]
dfs_in_pore = dfs.loc[dfs['realz'] < 0]

#write files
dfl_in_pore.to_csv('long_particles_in_pore.txt', sep=' ', header=None)
dfs_in_pore.to_csv('short_particles_in_pore.txt', sep=' ', header=None)

#3_contacts_in_pore
#load file
fread1 = open('long_particles_in_pore.txt', 'r')
longinitial = fread1.readlines()
fread2 = open('short_particles_in_pore.txt', 'r')
shortinitial = fread2.readlines()
# strip and split data
long = []
short = []
for i in longinitial:
    first = i.strip(',')
    second = first.split()
    long.append(second)
for i in shortinitial:
    first = i.strip(',')
    second = first.split()
    short.append(second)

#whether whole chain in pore
j = 0
kl = 0
ks = 0
lcontact = []
scontact = []
lchainnumber = int(long[0][1])
schainnumber = int(short[0][1])
longchainlength = 120 #should be changed
shortchainlength = 40 #should be changed
for i in range(0, len(long)-1,1):
    if int(long[i][1]) == int(lchainnumber):
        j = j + 1
    else:
        lcontact.append(j)
        j = 1
        lchainnumber = int(long[i][1])
j = 0
for i in range(0, len(short)-1,1):
    if int(short[i][1]) == int(schainnumber):
        j = j + 1
    else:
        scontact.append(j)
        j = 1
        schainnumber = int(short[i][1])

#write into a new file
fwrite1 = open('long_contacts_in_pore.txt', 'w')
fwrite2 = open('short_contacts_in_pore.txt', 'w')
for i in range(0, len(lcontact), 1):
    s = str(lcontact[i]).replace('[','').replace(']','')
    s = s.replace("'",'').replace(',','') +'\n'
    fwrite1.write(s)
for i in range(0, len(scontact), 1):
    a = str(scontact[i]).replace('[','').replace(']','')
    a = a.replace("'",'').replace(',','') +'\n'
    fwrite2.write(a)
fwrite1.close()
fwrite2.close()
