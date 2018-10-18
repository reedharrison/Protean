import os as os
import sys as sys
sys.path.append('../../')

import numpy as np
import mdtraj as md
from joblib import Parallel, delayed
import matplotlib.pyplot as plt

# from protean.optimize.SimulatedAnnealing import SimulatedAnnealing
from protean.editor.mutate import random_mutagenesis, refine_structures, generate_sites
# from protean.utilities.conversion import openmm2mdtraj, mdtraj2openmm

from test_sa import optimize
from simtk import unit

# def optimize(trj, moleculeSelections=['chainid 0 1', 'chainid 0', 'chainid 1']):
# 	topology, positions = mdtraj2openmm(trj)
# 	optimizer = SimulatedAnnealing(topology=topology, positions=positions, moleculeSelections=moleculeSelections)
# 	optimizer.run()
# 	return openmm2mdtraj(optimizer.getTopology(), optimizer.getPositions())

# 6.2 mins to generate and optimize 10 structures...

def main():
	pdbfile = 'PDL1-MbG9_400ps.pdb'
	nMutants = 10

	parent = md.load(pdbfile)

	# sites = generate_sites(parent, selection='chainid 0')
	query = parent.top.select('chainid 1')
	haystack = parent.top.select('chainid 0 and resSeq 74 to 93')
	indices = md.compute_neighbors(traj=parent, cutoff=0.5, query_indices=query, haystack_indices=haystack)
	sites = generate_sites(parent, indices=indices[0].tolist())

	with Parallel(n_jobs=4, verbose=50) as parallel:
		mdls = parallel(delayed(random_mutagenesis)(parent=parent, sites=sites) for x in range(nMutants))
		mdls = parallel(delayed(refine_structures)(trjs=mdl) for mdl in mdls)
		mdls = parallel(delayed(optimize)(trj=mdl) for mdl in mdls)
	return mdls

result = main()
mdls = [x[0] for x in result]
scores = [x[1].value_in_unit(unit.kilojoule_per_mole) for x in result]

ref = md.load('PDL1-MbG9_400ps.pdb')
ref_atoms = ref.top.select('name CA and chainid 0')

rmsds = []
for i, mdl in enumerate(mdls):
	mdl_atoms = mdl.top.select('name CA and chainid 0')
	mdl = mdl.superpose(ref, atom_indices=mdl_atoms, ref_atom_indices=ref_atoms)
	# val = md.rmsd(mdl, ref, atom_indices=mdl_atoms, ref_atom_indices=ref_atoms)

	xyz_ref = ref.atom_slice(ref_atoms)
	xyz_mdl = mdl.atom_slice(mdl_atoms)

	val = np.sqrt(np.square((xyz_mdl.xyz*10.) - (xyz_ref.xyz*10.)).sum() / float(xyz_ref.n_atoms * 3))

	rmsds.append(val)
	mdl.save('mdl_%d.pdb' % i)

fig, ax = plt.subplots(1,1)
ax.plot(rmsds, 'o')
ax.set_xlabel(r'Mutant Index')
ax.set_ylabel(r'RMSD ($\AA$)')
fig.savefig(r'Mutant_RMSD_from_Parent.pdf', dpi=600)

fig, ax = plt.subplots(1,1)
ax.plot(scores, rmsds, 'o')
ax.set_xlabel(r'Score ($kJ/mol$)')
ax.set_ylabel(r'RMSD ($\AA$)')
fig.savefig(r'Mutant_RMSD_vs_Score.pdf', dpi=600)

fig, ax = plt.subplots(1,1)
for i, (rmsd, score) in enumerate(zip(rmsds, scores)):
	ax.plot(score, rmsd, 'o', label=i)
ax.legend()
ax.set_xlabel(r'Score ($kJ/mol$)')
ax.set_ylabel(r'RMSD ($\AA$)')
fig.savefig(r'Mutant_RMSD_vs_Score.pdf', dpi=600)
