from optparse import OptionParser                                                                                                    
parser=OptionParser()
parser.add_option('-r', '--ref', dest='ref', help='reference cen180', default='')
parser.add_option("-a", "--acc", dest="acc", help="acc cen180", default="")
(options, args)=parser.parse_args()

import sys
from Bio import SeqIO

cens_acc={}
for cen in SeqIO.parse("%s"%options.acc, "fasta"):
    name, seq = str(cen.id) , hash(str(cen.seq))
    cens_acc[name] = seq


key_list = list(cens_acc.keys())
val_list = list(cens_acc.values())

out=open('%s.present'%options.acc, 'w')
out2=open('%s.absent'%options.acc, 'w')
cens_ref=set()
for record in SeqIO.parse('%s'%options.ref, "fasta"):
    seq_hash = hash(str(record.seq))
    cens_ref.add(hash(str(record.seq)))
    if seq_hash in cens_acc.values(): 
        position = val_list.index(seq_hash) 
        out.write(str(record.id) + '\t' + "1" + '\n')
    else: 
        out2.write(str(record.id) + '\t' + "0" + '\n')
