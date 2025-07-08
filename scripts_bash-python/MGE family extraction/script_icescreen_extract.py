import sys

infile = sys.argv[1]
# lista de familias
flist = []
# lista que vamos a contar
clist = []
#hacemos lista de contigs
cont_list = []
nl = 0
name_genom = str(sys.argv[1].split("/")[-1][:-16])
print (name_genom)

for line in open(infile):
	nl += 1
	if nl ==1:
		line = line.strip("\n")
		line = line.split("\t")
		fline = line[2] + "\t" + line[9]
	else:
		line = line.strip("\n")
		line = line.split("\t")
                # cogemos las familias del genoma (todas) en una lista
		clist.append(line[9])
                #cogemos los contigs en una lista
		if line[2] not in cont_list:
			cont_list.append(line[2])
                # cogemos cada tipo de familia del genoma (todos) en una lista
		if line[9] not in flist:
			flist.append(line[9])

# hacemos una lista final contando y asignando el conteo al tipo de familia
f_flist = []
for fam1 in clist:
	nf = clist.count(fam1)
	for fam2 in flist:
		if fam2 == fam1:
			entry = (fam2 + ":" + str(nf))
			if entry not in f_flist:
				f_flist.append(entry)

ofile = open("ices_family_count.tsv", "a")
# escribimos la linea del genoma
gline = (name_genom + "\t" + "\t".join(f_flist) + "\n")
ofile.write(gline)
#para cada contig

for cont in cont_list:
	clist = []
	flist = []
	for line in open(infile):
		line = line.strip("\n")
		line = line.split("\t")
		if line[2] == cont:
			clist.append(line[9])
			if line[9] not in flist:
				flist.append(line[9])
	f_flist = []
	for fam1 in clist:
		nf = clist.count(fam1)
		for fam2 in flist:
			if fam2 == fam1:
				entry = (fam2 + ":" + str(nf))
				if entry not in f_flist:
					f_flist.append(entry)

	ofile.write(cont + "\t" + "\t".join(f_flist) + "\n")
ofile.close()
