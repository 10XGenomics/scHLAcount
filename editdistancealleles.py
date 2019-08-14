import editdistance
import sys
from subprocess import Popen, PIPE

HLANUC="hla_nuc.fasta"

if len(sys.argv) != 3:
	print("Usage: python editdistancealleles.py A*01:01:01 A*03:02:01\nNote: if you do not include the maximum level of specificity (e.g. four fields) an ARBITRARY allele will be chosen matching the prefix given.")
	sys.exit()

a1 = sys.argv[1]
a2 = sys.argv[2]

a1p = Popen(["grep", "-F", "-m", "1", a1, HLANUC],stdout=PIPE,stderr=PIPE) 
(a1_stdout, _) = a1p.communicate()
a1_name = a1_stdout.split()[0][1:]
a1_allele = a1_stdout.split()[1]
a1p = Popen(["samtools", "faidx", HLANUC, a1_name],stdout=PIPE,stderr=PIPE) 
(a1_stdout, _) = a1p.communicate()
a1_seq = "".join([s.strip() for s in a1_stdout.split()[1:]])

a2p = Popen(["grep", "-F", "-m", "1", a2, HLANUC],stdout=PIPE,stderr=PIPE) 
(a2_stdout, _) = a2p.communicate()
a2_name = a2_stdout.split()[0][1:]
a2_allele = a2_stdout.split()[1]
a2p = Popen(["samtools", "faidx", HLANUC, a2_name],stdout=PIPE,stderr=PIPE) 
(a2_stdout, _) = a2p.communicate()
a2_seq =  "".join([s.strip() for s in a2_stdout.split()[1:]])

d = editdistance.eval(a1_seq,a2_seq)

print("CDS edit distance:", a1_allele, a2_allele, d)
