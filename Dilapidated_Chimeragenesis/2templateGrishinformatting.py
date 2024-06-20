#All formatting based on grishin format described in rosetta comparative modeling tutorial 2016
#Only intended to be used on clustalw alignments and with only 2 sequences
import numpy as np
import glob
file2 =open(glob.glob("*onSARS2.aln")[0],"r")
file =open(glob.glob("*onSARS.aln")[0],"r")
SARSAlignment=file.readlines()
SARS2Alignment=file2.readlines()
GrishinAlginment=np.empty((5,3), dtype=object)
GrishinAlginment2=np.empty((5,3), dtype=object)
Targetsequence=file.name.replace('onSARS.aln','')
Targetsequence2=file2.name.replace('onSARS2.aln','')
Templatesequence='SARS_INPUT'
Templatesequence2='SARS2_INPUT'
GrishinAlginment2[0,0]='##'
GrishinAlginment2[0,1]=Targetsequence2
GrishinAlginment2[0,2]=Templatesequence2 + '.pdb'
GrishinAlginment2[1,0]='#'
GrishinAlginment2[2,0]='scores from program : 0'
GrishinAlginment2[3,0]='0'
GrishinAlginment2[4,0]='0'
GrishinAlginment2[3,1]=''
GrishinAlginment2[4,1]=''
GrishinAlginment2[1,1]=''
GrishinAlginment2[1,2]=''
GrishinAlginment2[2,1]=''
GrishinAlginment2[2,2]=''
GrishinAlginment2[3,2]=''
GrishinAlginment2[4,2]=''
i=0
#4 lines separate the sequences in the clustalw alignment
for line in SARS2Alignment:
    i+=1
    if i%4==0:
        x=line.split()
        GrishinAlginment2[3,1]+=str(x[1])
    elif (i-1)%4==0 and i!=1:
        x = line.split()
        GrishinAlginment2[4, 1] += x[1]

GrishinAlginment[0,0]='##'
GrishinAlginment[0,1]=Targetsequence
GrishinAlginment[0,2]=Templatesequence+'.pdb'
GrishinAlginment[1,0]='#'
GrishinAlginment[2,0]='scores from program : 0'
GrishinAlginment[3,0]='0'
GrishinAlginment[4,0]='0'
GrishinAlginment[3,1]=''
GrishinAlginment[4,1]=''
GrishinAlginment[1,1]=''
GrishinAlginment[1,2]=''
GrishinAlginment[2,1]=''
GrishinAlginment[2,2]=''
GrishinAlginment[3,2]=''
GrishinAlginment[4,2]=''
i=0
#4 lines separate the sequences in the clustalw alignment
for line in SARSAlignment:
    i+=1
    if i%4==0:
        x=line.split()
        GrishinAlginment[3,1]+=str(x[1])
    elif (i-1)%4==0 and i!=1:
        x = line.split()
        GrishinAlginment[4, 1] += x[1]
Title=Targetsequence+'on'+Templatesequence+'.'+'grishin'
Title2=Targetsequence2+'on'+Templatesequence2+'.'+'grishin'
np.savetxt(Title,GrishinAlginment,fmt="%s",delimiter=" ")
np.savetxt(Title2,GrishinAlginment2,fmt="%s",delimiter=" ")
exit()