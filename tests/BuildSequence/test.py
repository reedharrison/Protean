import os as os
import sys as sys
sys.path.append('../../')

import mdtraj as md
from protean.model.homology import build_sequence
from protean.editor.mutate import refine_structures

# sequence is 
# sequence is from human HSF1 RD domain (>sp|Q00613|221-310) as defined by Uniprot
sequence = 'SMPKYSRQFSLEHVHGSGPYSAPSPAYSSSSLYAPDAVASSGPIISDITELAPASPMASPGGSIDERPLSSSPLVRVKEEPPSPPQSPRV'

trj = build_sequence(sequence)
trj = refine_structures(trj)
trj.save('HSF1_human_RD.pdb')
