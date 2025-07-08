# metemos la tabla con el conteo de familias para cada contig

import sys
infile = sys.argv[1]

#hacemos lista de familias
f_list = []
for line in open(infile):
	line = line.strip("\n")
	line = line.split("\t")
	line = sorted(line[1:])
	for element in line:
		fam = element.split(":")[0]
		if fam not in f_list:
			f_list.append(fam)
	contig = line[0]

f_list = sorted(f_list)
f_line = "contig" + "\t" + "\t".join(f_list) + "\n"

#cambiar nombre del file este era una prueba!
ofile = open("isescan_family.tsv","w")
ofile.write(f_line)
for line in open(infile):
	line = line.strip("\n")
	line = line.split("\t")
	nl = line[0] + "\t"
	line = sorted(line[1:], key=lambda x: x.split(':')[0])
	pmatch = 0
	for i in range(0,len(f_list)):
		for pair in line:
			if pair.split(":")[0] == f_list[i]:
				umatch = i
				esp = (umatch - pmatch)
				nl += (esp * ("Nan" + "\t"))
				nl += str(pair.split(":")[1] + "\t")
				pmatch = i + 1
			else :
				continue
		if (i == (len(f_list)-1)) and (pmatch != len(f_list)):
			esp = (i - pmatch)
			nl += (esp * ("Nan" + "\t"))
			nl += ("Nan")
	nl += "\n"
	ofile.write(nl)
ofile.close()
