import os as os
import sys as sys
sys.path.append('../../')

import pickle as p
import numpy as np
import mdtraj as md

from simtk import unit
from protean.evolve.evolve import Evolve

if __name__ is "__main__":
	pdbfile = 'c3d-s1p2_nopatches.pdb'

	# pdb = md.load(pdbfile)
	# query = pdb.top.select('chainid 1')
	# haystack = pdb.top.select('chainid 0 and resSeq 74 to 93')
	# indices = md.compute_neighbors(traj=pdb, cutoff=1., query_indices=query, haystack_indices=haystack)[0].tolist()

	optimizationOpts = {
					'temperatureRange': (350., 250.), 
					'temperatureSteps': (3, 50), 
					'simulationSteps': 100,
					'forcefield': None,
					'constructIntegrator': None, 
					'constructSystem': None, 
					'platform': None,
					'moleculeSelections':['chainid 0 1', 'chainid 0', 'chainid 1'], 
					'constSimulationTemp': None, 
					'constrainAtoms': None, 
					'restrainAtoms': None,
					'restraintConstant': 5.0*unit.kilojoules/(unit.angstrom**2)
	}

	evolver = Evolve(parent=pdbfile, nGenerations=5, nChildren=4, selection='chainid 0 and resSeq 1 to 9', nThreads=4,
		optimizationOpts=optimizationOpts)

	# evolver._generateChildren(0, verbose=50)
	# survivors = evolver._findSurvivors(0)
	# print('These children survive:')
	# print(survivors)

	evolver.run(verbose=50)
	p.dump(evolver, open('evolver_object.p', 'wb'))

# if __name__ is "__main__":
# 	evolver = main()
# 	scores = evolver.scores[:,:]
# 	p = boltzmann_p(scores)
# 	idx = np.argmax(p)
# 	idx = np.unravel_index(idx, scores.shape)
# 	print(idx)
