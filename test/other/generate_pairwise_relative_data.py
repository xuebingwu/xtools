
l2n = {'A':1.0, 'C':2.0, 'G':3.0, 'T':4.0}

f=open('')
for line in f:
	flds = line.strip().split()
	line2 = f.readline()
    flds2 = line2.strip().split()
	print flds[0]+"\t"+str(l2n[flds[0][29]]/l2n[flds2[0][29]])
    print flds[0]+"\t1.0"

f.close()


