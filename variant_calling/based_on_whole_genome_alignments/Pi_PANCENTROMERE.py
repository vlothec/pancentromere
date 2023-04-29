from optparse import OptionParser
parser = OptionParser()
parser.add_option("-c", "--csv", dest="csv", help="combined csv file", default="")
parser.add_option("-d", "--dir", dest="dir", help="directory of combined vcf file", default="")

(options, args) = parser.parse_args()

out=open('%s/%s.distance'%(options.dir, options.csv), 'w')

from itertools import combinations

for line in open('%s/%s'%(options.dir,options.csv),'r'):
    line=line.strip().split(',')
    if 'chr' in line:
        header=line
        acc=line[2:]  
        comb=list(combinations(acc,2))
        #print(comb)
        #break
        #print(len(comb))
        #break
        totalL,dist=[0]*len(comb),[0]*len(comb) 
    else:
        gen=line[2:]	
        gen=list(filter(lambda a: a != 'NA', gen))
        nonV,V=0,0	
        if len(set(gen))==1:
            nonV=1
        if len(list(set(gen)))>1:
            V=1	
        count=0
        for i,j in list(combinations(acc,2)):
            ind_i, ind_j=acc.index(i),acc.index(j)
            gen_i,gen_j=gen[ind_i],gen[ind_j]
            if gen_i==gen_j:
                totalL[count]+=1
            else: 
                dist[count]+=abs(float(gen_i)-float(gen_j))
                #if V==1:
                #    dist[count]+=abs(float(gen_i)-float(gen_j))
            count+=1
out.write('acc1\t%acc2\taln\ttotal_dist\t\n')
for i in range(0, len(comb)):
    out.write('%s\t%s\t%s\t%s\n'%(list(comb[i])[0], list(comb[i])[1], totalL[i], dist[i]))
out.close()

