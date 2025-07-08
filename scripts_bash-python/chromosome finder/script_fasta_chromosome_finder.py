#cogemos un fasta metemos contig y  su secuencia en un diccionario.
# Al contig mas largo le llamamos cromosoma
# hacemos un tsvo que tenga como fila el nombre del genoma y debajo sus contigs con el cromosoma indicado.


import sys
infile = sys.argv[1]

genome_name = str(infile.split("/")[-1][:-4])
print (genome_name)

dicont = {}
n_contig = 0
seq = ''
for line in open (infile):
	if line[0] == ">":
		n_contig += 1
		line = line.strip("\n")
		contig_name = line.strip(">")
		print (contig_name)
		seq = ''
		contig_n = "seq_" + str(n_contig)
	else:
		line = line.strip("\n") 
		seq += line
		dicont[contig_name] = seq

candidate_l = 0
chrom_l = 0

for key in dicont.keys():
	candidate_l = len (dicont[key])
	print (candidate_l)
	if candidate_l > chrom_l:
		chrom_l = candidate_l
		chrom = key

new_name = str(chrom) + "_chromosome"
dicont[new_name] = dicont[chrom]
del dicont[chrom]

print (new_name)

ofile = open("accesion_contig_chromosome_list.txt","a")
ofile.write (genome_name + "\n")

for key in sorted(dicont.keys()):
	ofile.write (key + "\n")
ofile.close()
