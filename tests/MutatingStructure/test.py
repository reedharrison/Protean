"""
Test Manual Mutations
=====================
"""

import os as os
import sys as sys
sys.path.append('../../')

import mdtraj as md

from protean.editor.mutate import mutate_structure, refine_structures, random_mutagenesis

pdbfile = 'PDL1-MbG9_400ps.pdb'
trj = md.load(pdbfile)

mutations = [{'chainid': 0, 'position': 78, 'type': 'A'}]

mdl = mutate_structure(parent=trj, mutations=mutations)
mdl.save('mutant_raw.pdb')

mdl = refine_structures(mdl, minimize=True)
mdl.save('mutant_refined.pdb')


"""
Test Random Mutagenesis
=======================
"""
import os as os
import sys as sys
sys.path.append('../../')

import mdtraj as md
from random import random
from protean.editor.mutate import mutate_structure, refine_structures, random_mutagenesis

pdbfile = 'PDL1-MbG9_400ps.pdb'
trj = md.load(pdbfile)

sites = []
prob = 0.1
for i, chain in enumerate(trj.top.chains):
	for j, residue in enumerate(chain.residues):
		if random() <= prob:
			sites.append(tuple([i, j]))

mdl = random_mutagenesis(parent=trj, sites=sites)
mdl.save('random_mutant_raw.pdb')

mdl = refine_structures(mdl, minimize=True)
mdl.save('random_mutant_refined.pdb')
