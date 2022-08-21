from optparse import OptionParser
parser = OptionParser()
parser.add_option("-v", "--vcf", dest="vcf", help="combined vcf file", default="")
parser.add_option("-d", "--dir", dest="dir", help="directory of combined vcf file", default="")
(options, args) = parser.parse_args()

out=open('%s/%s.csv'%(options.dir, options.vcf), 'w')
for line in open('%s/%s'%(options.dir, options.vcf), 'r'): 
	if '##' in line:
		pass
	elif '#CHR' in line:
		names=line.strip().split()[9:]
		out.write('chr,pos,%s\n'%(','.join(names)))
	elif '#' not in line:
		gen=line.strip().split()[9:]
		gen=[x.split(':')[0] for x in gen]
		gen_new=[]
		for g in gen:
			if g=='./.' or g =='./././.':
				gen_new.append('NA')
			else:
				g=[int(x) for x in g.split('/')]
				gen_new.append(str(sum(g)/float(len(g))))
		test=filter(lambda a: a != 'NA', gen_new)
		if len(set(test))>0:
			gen_new=','.join(gen_new)
			out.write('%s,%s,%s\n'%(line.strip().split()[0], line.strip().split()[1], gen_new))
	else:
		pass
