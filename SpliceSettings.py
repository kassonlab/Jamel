FastaSplicingInputs={"Fasta1":'FullSARS2.fasta','Fasta2':'AvianD274.fasta','Boundary1':630,'Boundary2':635,'NumofShifts':127,'ShiftLength':5,'MobileBoundary':'Left'}
from Chimeragenesis import ShiftedSplice

ShiftedSplice.SplicingFastas(**FastaSplicingInputs)