from directional_bins import directional_consolidation

f1='test/data/complex/r1_out.fastq'
f2='test/data/complex/r2_out.fastq'
o1='test/data/complex/r1_cons.fastq'
o2='test/data/complex/r2_cons.fastq'

res = {}
for i in range(20):
    res[i + 1] = directional_consolidation(f1, f2, o1, o2, 0, .6, i+1)
print res

