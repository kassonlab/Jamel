from numpy import empty,savetxt
from Analysis import relative_stability,overall_confidence,averaging_multimer_plddt

protein_list=[line.split() for line in open('/gpfs/gpfs0/scratch/jws6pq/BridCMfiles/Flu_CovidList', 'r').readlines()]
PDB=[x +'.pdb' for x in protein_list]
plddt_files=[x + '.plddt' for x in protein_list]

averaged_plddt=list(map(averaging_multimer_plddt, plddt_files))

splice_length_one=200
boundary_one=[x for x in range(0, 489 - splice_length_one, 10)]
protein_one=['HA.fasta' for x in boundary_one]
boundary_two=[x + splice_length_one for x in boundary_one]

splice_length_two=200
boundary_three=[x for x in range(540, 973 - splice_length_two, 20)]
protein_two=['SARS2.fasta' for x in boundary_three]
boundary_four=[x + splice_length_one for x in boundary_three]

relative_stability_list=[]
HA_relative=[]
SARS_relative=[]
k=0
for i in range(len(boundary_one)):
    for j in range(len(boundary_three)):
        protein_one_rs = relative_stability('AvgHA.plddt', averaged_plddt[k], [0, 0 + splice_length_one], [boundary_one[i], boundary_two[i]])
        protein_two_rs = relative_stability('Avg3merSARS2.plddt', averaged_plddt[k], [0 + splice_length_one, None], [boundary_three[j], boundary_four[j]])
        relative_stability_list.append((protein_one_rs[0] + protein_two_rs[0]) / (protein_one_rs[1] + protein_two_rs[1]))
        HA_relative.append(protein_one_rs[0] / protein_one_rs[1])
        SARS_relative.append(protein_two_rs[0] / protein_two_rs[1])
        k+=1
overall_stability=list(map(overall_confidence, plddt_files))

data_chart=empty((len(protein_list) + 1, 5), dtype=object)
data_chart[0, 0], data_chart[1:, 0]= 'Protein', protein_list
data_chart[0, 1], data_chart[1:, 1]= 'Average Stability Difference', relative_stability_list
data_chart[0, 2], data_chart[1:, 2]= 'HA Difference', HA_relative
data_chart[0, 3], data_chart[1:, 3]= 'SARS Difference', SARS_relative
data_chart[0, 4], data_chart[1:, 4]= 'Overall chimera plddt', overall_stability


savetxt('/gpfs/gpfs0/scratch/jws6pq/BridCMfiles/120522HACovidAnalysis.tsv', data_chart, fmt="%s,%s,%s,%s,%s", delimiter="")
