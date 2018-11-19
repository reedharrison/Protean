from __future__ import print_function, division

import numpy as np
import mdtraj as md
import simtk.unit as unit
from random import choice, random

from protean.evolve.defaults import mutant_filename_constructor
from protean.optimize.SimulatedAnnealing import SimulatedAnnealing
from protean.editor.mutate import random_mutagenesis, random_recombination, refine_structures, get_sequence
from protean.utilities.conversion import openmm2mdtraj, mdtraj2openmm
from protean.model.homology import homology_model


def _variationKernel(filenames, generation_index, child_index, sites, #seq_db, score_db,
	firstParent=False, percentRecombination=-1.0,
	minimization=True, tolerance=5.*unit.kilojoule_per_mole, maxIterations=1000,
	mutationOpts=None, refinementOpts=None, optimizationOpts=None, 
	retain_models=True, workdir=None, maxAttempts=5):

	assert (filenames is not None) or (workdir is not None)
	if (firstParent is not False) and (generation_index != 0):
		filenames = filenames + [firstParent]

	# Prepare
	energy = np.nan
	sequence = None
	if mutationOpts is None:
		mutationOpts = {'degree':1, 'library': None}
	if refinementOpts is None:
		refinementOpts = {
			'variants':None, 
			'forcefield':None, 
			'platform':None, 
			'pH':7.0, 
			'minimize':True,
			'maxIterations':1000,
			'energyTol':5.0*unit.kilojoule_per_mole
		}
	if optimizationOpts is None:
		optimizationOpts = {
			'temperatureRange': (350., 250.), 
			'temperatureSteps': (11, 50), 
			'simulationSteps': 100,
			'forcefield': None,
			'constructIntegrator': None, 
			'constructSystem': None, 
			'platform': None,
			'moleculeSelections':['protein'], 
			'constSimulationTemp': None, 
			'constrainAtoms': None, 
			'restrainAtoms': None,
			'restraintConstant': 5.0*unit.kilojoules/(unit.angstrom**2)
		}

	# Find starting sequence
	if (child_index == 0) and (firstParent is not False):
		sequence = get_sequence(md.load(firstParent))
	elif (random() <= percentRecombination) and (len(filenames) >= 2):
		sequence = random_recombination(parents=[get_sequence(md.load(x)) for x in filenames])
	else:
		sequence = get_sequence(md.load(choice(filenames)))

	if (firstParent is not False) or (child_index != 0):
		sequence = random_mutagenesis(parent=sequence, sites=sites, **mutationOpts)

	mdl = _homologyKernel(sequence=sequence, templates=filenames)
	mdl = _refinementKernel(mdl, **refinementOpts)
	mdl, score = _optimizationKernel(mdl, minimization=minimization, tolerance=tolerance, 
		maxIterations=maxIterations, **optimizationOpts)

	if retain_models:
		if workdir is not None:
			outpath = workdir
		else:
			outpath = os.path.dirname(filename)
		pdbout = mutant_filename_constructor(outpath, generation_index, child_index)
		mdl.center_coordinates().save(pdbout)

	energy = score.value_in_unit(unit.kilojoule_per_mole)
	sequence = get_sequence(mdl)

	return (energy, sequence)

# def _variationKernel(filenames, generation_index, child_index, sites, #seq_db, score_db, # need to replace method to support recombination
# 	firstParent=False, percentRecombination=-1.0,
# 	minimization=True, tolerance=5.*unit.kilojoule_per_mole, maxIterations=1000,
# 	mutationOpts=None, refinementOpts=None, optimizationOpts=None, 
# 	retain_models=True, workdir=None, maxAttempts=5):

# 	assert (filenames is not None) or (workdir is not None)
# 	if (firstParent is not False) and (generation_index != 0):
# 		filenames = filenames + [firstParent]

# 	# Prepare
# 	energy = np.nan
# 	sequence = None
# 	if mutationOpts is None:
# 		mutationOpts = {'degree':1, 'library': None}
# 	if refinementOpts is None:
# 		refinementOpts = {
# 			'variants':None, 
# 			'forcefield':None, 
# 			'platform':None, 
# 			'pH':7.0, 
# 			'minimize':True,
# 			'maxIterations':1000,
# 			'energyTol':5.0*unit.kilojoule_per_mole
# 		}
# 	if optimizationOpts is None:
# 		optimizationOpts = {
# 			'temperatureRange': (350., 250.), 
# 			'temperatureSteps': (11, 50), 
# 			'simulationSteps': 100,
# 			'forcefield': None,
# 			'constructIntegrator': None, 
# 			'constructSystem': None, 
# 			'platform': None,
# 			'moleculeSelections':['protein'], 
# 			'constSimulationTemp': None, 
# 			'constrainAtoms': None, 
# 			'restrainAtoms': None,
# 			'restraintConstant': 5.0*unit.kilojoules/(unit.angstrom**2)
# 		}

# 	# Try to generate model
# 	counter = 0
# 	while counter < maxAttempts:

# 		# Find starting sequence
# 		if (child_index == 0) and (firstParent is not False):
# 			sequence = get_sequence(md.load(firstParent))
# 		elif (random() <= percentRecombination) and (len(filenames) >= 2):
# 			sequence = random_recombination(parents=[get_sequence(md.load(x)) for x in filenames])
# 		else:
# 			sequence = get_sequence(md.load(choice(filenames)))

# 		if (firstParent is not False) or (child_index != 0):
# 			sequence = random_mutagenesis(parent=sequence, sites=sites, **mutationOpts)

# 		mdl = _homologyKernel(sequence=sequence, templates=filenames)

# 		try:
# 			mdl = _refinementKernel(mdl, **refinementOpts)
# 			mdl, score = _optimizationKernel(mdl, minimization=minimization, tolerance=tolerance, 
# 				maxIterations=maxIterations, **optimizationOpts)

# 			if retain_models:
# 				if workdir is not None:
# 					outpath = workdir
# 				else:
# 					outpath = os.path.dirname(filename)
# 				pdbout = mutant_filename_constructor(outpath, generation_index, child_index)
# 				mdl.center_coordinates().save(pdbout)

# 			# seq_db[generation_index, child_index] = get_sequence(mdl)
# 			# score_db[generation_index, child_index] = score.value_in_unit(unit.kilojoule_per_mole)

# 			energy = score.value_in_unit(unit.kilojoule_per_mole)
# 			sequence = get_sequence(mdl)

# 			counter = maxAttempts

# 		except:
# 			print('\nERROR: failed to generate child %d from generation %d ... attempt %d\n' % (child_index, generation_index, counter))
# 			# seq_db[generation_index, child_index] = None
# 			# score_db[generation_index, child_index] = np.nan

# 		counter += 1

# 	return (energy, sequence)

# def _mutagenesisKernel(trj, sites, **kwargs):
# 	"""
# 	inputs
# 	======
# 	trj: mdtraj.Trajectory
# 		Structure to be mutated. Must be object and not filename in path.
# 	sites: list
# 		List of mutatable positions of the form (<chain index>, <residue position index>).
# 		All indices start at 0 for each chain. Thus, neither chain index or position
# 		index will be unique. Only the tuple of chain index and position index will be
# 		unique.
# 	**kwargs:
# 		Kewword arguments suitable for protean.mutate.mutate.random_mutagenesis

# 	outputs
# 	=======
# 	out: mdtraj.Trajectory
# 		Mutated structure

# 	"""
# 	newtrj = random_mutagenesis(trj, sites, **kwargs)
# 	return newtrj.center_coordinates()

def _homologyKernel(sequence, templates, **kwargs):
	if isinstance(templates, md.Trajectory):
		templates = [templates]
	elif isinstance(templates, str):
		templates = [md.load(templates)]
	elif isinstance(templates, list):
		if all([isinstance(x, str) for x in templates]):
			templates = [md.load(x) for x in templates]

	mdl = homology_model(sequence, templates, **kwargs)
	# mdl.top.create_disulfide_bonds(mdl.openmm_positions(-1)) # 10/11/2018, try to update bond block for DISU to prevent topology errors with CYS
	return mdl.center_coordinates()

def _refinementKernel(trj, **kwargs):
	"""
	inputs
	======
	trj: mdtraj.Trajectory
		Input structure to refine
	**kwargs:
		Keyword arguments suitable for protean.model.addH.protonate_protein

	outputs
	=======
	out: mdtraj.Trajectory
		Refined structure
	"""
	newtrj = refine_structures(trj, **kwargs)
	return newtrj.center_coordinates()

def _optimizationKernel(trj, minimization=True, tolerance=5.*unit.kilojoule_per_mole, 
	maxIterations=1000, **kwargs):
	"""
	inputs
	======
	trj: mdtraj.Trajectory
		Input structural file to be optimized
	**kwargs: dict
		Extra keyword arguments suitable for SimulatedAnnealing class

	outputs
	=======
	out: tuple
		Tuple containting new trajectory (mdtraj.Trajectory) and associated optimization 
		score (simtk.unit.Quantity) based on the potential energy of the structure. 
	"""
	topology, positions = mdtraj2openmm(trj)
	optimizer = SimulatedAnnealing(topology=topology, positions=positions, **kwargs)
	if minimization:
		optimizer.minimize(maxIterations=maxIterations, tolerance=tolerance)
	optimizer.run()
	newtrj = openmm2mdtraj(optimizer.getTopology(), optimizer.getPositions())
	score = optimizer.getScore()
	return (newtrj.center_coordinates(), score)

def _retrieveStructure(workdir, generation_index, child_index):
	filename = mutant_filename_constructor(workdir, generation_index, child_index)
	return md.load(filename)

def _generateStructure(sequence, templates, homologyOpts=None,
	minimization=True, tolerance=5.*unit.kilojoule_per_mole, maxIterations=1000,
	mutationOpts=None, refinementOpts=None, optimizationOpts=None):

	if homologyOpts is None:
		homologyOpts = {'env':None, 'n_models':1, 'return_all_models':False}
	if refinementOpts is None:
		refinementOpts = {'variants':None, 'forcefield':None, 'platform':None, 'pH':7.0, 'minimize':False}
	if optimizationOpts is None:
		optimizationOpts = {
			'temperatureRange': (350., 250.), 
			'temperatureSteps': (11, 50), 
			'simulationSteps': 100,
			'forcefield': None,
			'constructIntegrator': None, 
			'constructSystem': None, 
			'platform': None,
			'moleculeSelections':['protein'], 
			'constSimulationTemp': None, 
			'constrainAtoms': None, 
			'restrainAtoms': None,
			'restraintConstant': 5.0*unit.kilojoules/(unit.angstrom**2)
		}

	mdl = _homologyKernel(sequence, templates, **homologyOpts)
	mdl = _refinementKernel(mdl, **refinementOpts)
	mdl, score = _optimizationKernel(mdl, minimization=minimization, tolerance=tolerance, 
		maxIterations=maxIterations, **optimizationOpts)
	return mdl.center_coordinates()
