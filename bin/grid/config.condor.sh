


# Requirements = ( (machine == "macfrank.rice.edu") )
Rank = ( (machine == "macfrank.rice.edu")*3 + (machine == "star3.local")*2 + (machine == "star4.local")*2 )

InitialDir = /home/jdb12/work/dimuonAna/FemtoDstPairAna/bin/
Executable = /home/jdb12/work/dimuonAna/FemtoDstPairAna/bin/pairAna.app

GetEnv     = True

Queue 96