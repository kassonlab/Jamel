
FastaExchangeInputs={"Fastafile1":'FullSARS2.fasta','Fastafile2':'AvianD274.fasta','Boundary1':630,'Boundary2':635,"Fastalist":''}
file3 = open(FastaExchangeInputs['Fastalist'], "r")
list=file3.readlines()
import RBDExchanger
AlphaFoldEntry=''
for line in list:
   FastaExchangeInputs['Fastafile2']=line.rstrip()+'.fasta'
   AlphaFoldEntry+=(RBDExchanger.RBDExchange(FastaExchangeInputs['Fastafile1'],FastaExchangeInputs['Fastafile2'],FastaExchangeInputs['Boundary1'],FastaExchangeInputs['Boundary2'])) + ','
print(AlphaFoldEntry)
