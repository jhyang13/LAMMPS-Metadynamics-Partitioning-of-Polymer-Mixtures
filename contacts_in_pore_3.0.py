#一次处理多个文档
#包含所有的configuration

import pandas as pd

#1_rearrange
#load file
EPS = ['0.20', '0.40', '0.60', '0.80', '1.00', '1.20', '1.40', '1.60', '1.80', '2.00']
REP = ['01', '02', '03']

for z in EPS:
    for i in REP:
        fread1 = open('D:\\all\\projects\\python_research_code\\py_count_contacts_in_pore\\20221220_NA20_PHIA0.425_NB5_PHIB0.425_R2_EPS' + z + '_REP' + i + '\\long.lammpstrj', 'r')
        longread1 = fread1.readlines()
        fread2 = open('D:\\all\\projects\\python_research_code\\py_count_contacts_in_pore\\20221220_NA20_PHIA0.425_NB5_PHIB0.425_R2_EPS' + z + '_REP' + i + '\\short.lammpstrj', 'r')
        shortread1 = fread2.readlines()
        # strip and split data
        longread2 = []
        shortread2 = []
        for j in longread1:
            first = j.strip('/n')
            second = first.split()
            longread2.append(second)
        for j in shortread1:
            first = j.strip('/n')
            second = first.split()
            shortread2.append(second)
        #delete title, 8 means the length of one line in dump file
        longdump = []
        shortdump = []
        for j in range(0, len(longread2), 1):
            if len(longread2[j]) == 8:
                longdump.append(longread2[j])
        for j in range(0, len(shortread2), 1):
            if len(shortread2[j]) == 8:
                shortdump.append(shortread2[j])
        # add time step
        # parameters should be changed according to the dump files
        timesteplong = []
        timestepshort = []
        idumpstep = 101000000  # timestep the dump start
        dumpsegment = 100000  # dump every time step
        totalparticleslong = 6640  # atoms number
        totalparticlesshort = 6640  # atoms number
        for j in range(1, len(longdump) // totalparticleslong + 1, 1):
            for k in range(1, totalparticleslong + 1, 1):
                timesteplong.append(idumpstep)
            idumpstep = idumpstep + dumpsegment
        idumpstep = 101000000  # reset the initial timestep
        for j in range(1, len(shortdump) // totalparticlesshort + 1, 1):
            for k in range(1, totalparticlesshort + 1, 1):
                timestepshort.append(idumpstep)
            idumpstep = idumpstep + dumpsegment
        #add chain number
        longchainnumber = []
        shortchainnumber = []
        longchainlength = 20 #chain length
        shortchainlength = 5 #chain length
        totalparticleslong = 6640 #atoms number
        totalparticlesshort = 6640  #atoms number
        for j in range(1, len(longdump)//totalparticleslong+1, 1):
            p = 0
            for k in range(1, totalparticleslong//longchainlength+1, 1):
                p = p + 1
                for h in range(1, longchainlength+1, 1):
                    longchainnumber.append(p)
        for j in range(1, len(shortdump)//totalparticlesshort+1, 1):
            p = 0
            for k in range(1, totalparticlesshort//shortchainlength+1, 1):
                p = p + 1
                for h in range(1, shortchainlength+1, 1):
                    shortchainnumber.append(p)
        #write into a new file
        fwrite1 = open('long_rearrange.txt', 'w')
        fwrite2 = open('short_rearrange.txt', 'w')
        a = 'timestep chainN ID TYPE x y z ix iy iz ' + '\n'
        d = 'timestep chainN ID TYPE x y z ix iy iz ' + '\n'
        fwrite1.write(a)
        fwrite2.write(d)
        for j in range(0, len(longchainnumber), 1):
            a = str(timesteplong[j]).replace('[', '').replace(']', '')
            a = a.replace("'", '').replace(',', '') + ' '
            b = str(longchainnumber[j]).replace('[', '').replace(']', '')
            b = b.replace("'", '').replace(',', '') + ' '
            c = str(longdump[j]).replace('[','').replace(']','')
            c = c.replace("'",'').replace(',','') + '\n'
            fwrite1.write(a)
            fwrite1.write(b)
            fwrite1.write(c)
        for j in range(0, len(shortchainnumber), 1):
            d = str(timestepshort[j]).replace('[', '').replace(']', '')
            d = d.replace("'", '').replace(',', '') + ' '
            e = str(shortchainnumber[j]).replace('[', '').replace(']', '')
            e = e.replace("'", '').replace(',', '') + ' '
            f = str(shortdump[j]).replace('[','').replace(']','')
            f = f.replace("'",'').replace(',','') + '\n'
            fwrite2.write(d)
            fwrite2.write(e)
            fwrite2.write(f)
        fwrite1.close()
        fwrite2.close()

        #2_pandas_particles_in_pore
        #load files
        dfl = pd.read_table('long_rearrange.txt', delim_whitespace=True)
        dfs = pd.read_table('short_rearrange.txt', delim_whitespace=True)
        #calculate real coordination
        xbox = 25 #box size
        ybox = 25 #box size
        zbox = 36 #box size
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
        for j in longinitial:
            first = j.strip(',')
            second = first.split()
            long.append(second)
        for j in shortinitial:
            first = j.strip(',')
            second = first.split()
            short.append(second)
        #whether whole chain in pore
        k = 0
        lcontact = []
        ltimestep = []
        scontact = []
        stimestep = []
        lchainnumber = int(long[0][2])
        schainnumber = int(short[0][2])
        for j in range(0, len(long)-1,1):
            if int(long[j][2]) == int(lchainnumber):
                k = k + 1
            else:
                lcontact.append(k)
                ltimestep.append(int(long[j-1][1]))
                k = 1
                lchainnumber = int(long[j][2])
        k = 0
        for j in range(0, len(short)-1,1):
            if int(short[j][2]) == int(schainnumber):
                k = k + 1
            else:
                scontact.append(k)
                stimestep.append(int(short[j-1][1]))
                k = 1
                schainnumber = int(short[j][2])
        #write into a new file
        fwrite1 = open('long_contacts_in_pore_EPS' + z + '_' + i + '.txt', 'w')
        fwrite2 = open('short_contacts_in_pore_EPS' + z + '_' + i + '.txt', 'w')
        for j in range(0, len(lcontact), 1):
            a = str(ltimestep[j]).replace('[','').replace(']','')
            a = a.replace("'",'').replace(',','') +' '
            b = str(lcontact[j]).replace('[','').replace(']','')
            b = b.replace("'",'').replace(',','') +'\n'
            fwrite1.write(a)
            fwrite1.write(b)
        for j in range(0, len(scontact), 1):
            a = str(stimestep[j]).replace('[','').replace(']','')
            a = a.replace("'",'').replace(',','') +' '
            b = str(scontact[j]).replace('[','').replace(']','')
            b = b.replace("'",'').replace(',','') +'\n'
            fwrite2.write(a)
            fwrite2.write(b)
        fwrite1.close()
        fwrite2.close()