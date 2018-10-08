from __future__ import print_function, division

import os as os
import numpy as np
import zarr as zarr
import mdtraj as md
import simtk.unit as unit

from joblib import Parallel, delayed

# from protean.utilities.database import Database, SafeHDFStore
# from protean.optimize.SimulatedAnnealing import SimulatedAnnealing
# from protean.editor.mutate import random_mutagenesis, refine_structures, generate_sites
from protean.evolve.defaults import mutant_filename_constructor
from protean.editor.mutate import generate_sites
from protean.evolve.parallel import _generateStructure, _variationKernel
from protean.evolve.fitness import boltzmann_p


class Evolve:
	def __init__(self, parent, workdir=None, nGenerations=3, nChildren=100,
		mutationDegree=1, survivalCutoff=0.2, nThreads=-2, library=None,
		sites=None, selection='protein', atom_indices=None, retain_models=True,
		refinementOpts=None, optimizationOpts=None):
		"""
		Import parent structure
		=======================
		"""
		assert isinstance(parent, str) or isinstance(parent, md.Trajectory)
		if isinstance(parent, str):
			pdb = md.load(parent)
		elif isinstance(parent, md.Trajectory):
			pdb = parent

		"""
		Configure method
		================
		"""
		self._nGenerations = int(nGenerations)
		self._nChildren = int(nChildren)
		self._mutationDegree = int(mutationDegree)
		self._survivalCutoff = float(survivalCutoff)
		self._library = library
		self._nThreads = nThreads
		self._retain_models = retain_models

		if atom_indices is not None:
			self._sites = generate_sites(pdb, indices=atom_indices)
		else:
			self._sites = generate_sites(pdb, selection=selection)

		if refinementOpts is None:
			self._refinementOpts = {'variants':None, 'forcefield':None, 'platform':None, 'pH':7.0, 'minimize':False}
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
		self._optimizationOpts = optimizationOpts

		"""
		Initialize database
		===================
		For simplicity, I will skip this for now and just write out all structural files
			to  a working directory...
		"""
		# self.db = Database(database)
		# self._initializeDB()
		if workdir is None:
			workdir = os.path.join(os.getcwd(), 'protean_evolution')
		self._initialize_workdir(workdir)
		self._workdir = workdir

		scoreSynchronizer = zarr.ProcessSynchronizer(os.path.join(self._workdir, 'scores.sync'))

		# Persistence mode: ‘r’ means read only (must exist); ‘r+’ means read/write (must exist); 
		#	‘a’ means read/write (create if doesn’t exist); ‘w’ means create (overwrite if exists); 
		#	‘w-‘ means create (fail if exists).
		self.scores = zarr.open(os.path.join(self._workdir, 'scores.zarr'),
			mode='w', 
			shape=(self._nGenerations, self._nChildren),
			chunks=(1, self._nChildren),
			dtype=float,
			synchronizer=scoreSynchronizer
			)

		seqSynchronizer = zarr.ProcessSynchronizer(os.path.join(self._workdir, 'sequences.sync'))
		self.sequences = zarr.open(os.path.join(self._workdir, 'sequences.zarr'),
			mode='w',
			shape=(self._nGenerations, self._nChildren),
			chunks=(1, self._nChildren),
			dtype=str,
			synchronizer=seqSynchronizer
			)

		self.parent = os.path.join(self._workdir, 'parent.pdb')
		pdb.save(self.parent)

		self.survivors = [self.parent]

		"""
		Import parent structure
		=======================
		"""
		# inital_structures = []
		# if isinstance(parent, list):
		# 	for record in parent:
		# 		assert isinstance(record, str) or isinstance(record, md.Trajectory)
		# 		if isinstance(record, str):
		# 			pdb = md.load(record)
		# 		elif isinstance(record, md.Trajectory):
		# 			pdb = record
		# 		initial_structures.append(pdb)
		# else:
		# 	assert isinstance(parent, str) or isinstance(parent, md.Trajectory)
		# 	if isinstance(parent, str):
		# 		pdb = md.load(parent)
		# 	elif isinstance(parent, md.Trajectory):
		# 		pdb = parent
		# 	initial_structures.append(pdb)

		# self.parents = []
		# for i, structure in enumerate(initial_structures):
		# 	self.parents.append(os.path.join(self._workdir, 'parent_%d.pdb' % i))
		# 	structure.save(parent[-1])

		"""
		Add hidden variable to keep track if method has been run
		========================================================
		"""
		self._notRun = True

		return

	def _initialize_workdir(self, directory):
		if not os.path.exists(directory):
			os.makedirs(directory)
		return

	def _saveStructure(self, trj, generation_index, child_index):
		filename = mutant_filename_constructor(self._workdir, generation_index, child_index)
		trj.save(filename)
		return

	def _deleteStructure(self, generation_index, child_index):
		filename = mutant_filename_constructor(self._workdir, generation_index, child_index)
		if os.path.isfile(filename):
			os.remove(filename)
		else:
			print('Error: %s file not found!' % filename)
		return

	def _generateStructures(self, indices, templates=None, verbose=0, **kwargs):
		# hidden method to generate structure from sequence.
		# useful if structural file gets deleted. saves to work
		# directory for persistence.
		#
		# indicies is a tuple of the form (gen idx, child idx)
		sequences = [self.sequences[g, c] for g, c in indices]
		if templates is None:
			templates = self.parent

		trj = None
		with Parallel(n_jobs=self._nThreads, prefer='processes', verbose=verbose) as parallel:
			trj = parallel(delayed(_generateStructure)(sequence=x, templates=templates, **kwargs) for x in sequences)
		return trj

	def _deleteStructures(self, indices):
		# delete PDB files in work directory if they are not needed (AKA not survivors...)
		#
		# indices is a tuple of the form (gen idx, child idx)
		for g, c in indices:
			self._deleteStructure(generation_index=g, child_index=c)
		return

	def getStructure(self, indices, **kwargs):
		# retrieve structure from work directory, if possible.
		# if structure does not exist, generate a model.
		if self._notRun:
			print('Please call evolve.run() before attempting to retrieve structures!')
		else:
			models = {}
			keys_remodel = []
			trj_remodel = [] 
			for g, c in indices:
				filename = mutant_filename_constructor(self.workdir, generation_index=g, child_index=c)
				if os.path.isfile(filename):
					models[(g, c)] = md.load(filename)
				else:
					keys_remodel.append((g, c))
			trj_remodel = _generateStructures(indices=keys_remodel, templates=self.parent, **kwargs)
			for key, trj in zip(keys_remodel, trj_remodel):
				models.append(key, trj)
		return [models[x] for x in indices]

	def getSequence(self, indices):
		# retrieve sequence
		if (not all([(isinstance(x, list) or isinstance(x, tuple)) and (len(x)==2) for x in indices])) and (len(indices) == 2):
			indices = [indices]
		assert all([(isinstance(x, list) or isinstance(x, tuple)) and (len(x)==2) for x in indices])
		return [self.sequences[g, c] for g, c in indices]

	def _recombineSurvivors(self, generation_index, children_indices):
		# recombine child sequences to generate new children for the subsequent generation
		return

	def _generateChildren(self, generation_index, verbose=0):
		# parallel generation of children within a single generation
		# call recombine survivors here since natural variation should not
		# replace structures from recombination...
		with Parallel(n_jobs=self._nThreads, verbose=verbose) as parallel:
			parallel(delayed(_variationKernel)
				(
					filenames=self._findSurvivors(generation_index-1),
					workdir=self._workdir,
					generation_index=generation_index, 
					child_index=i, 
					sites=self._sites, 
					seq_db=self.sequences, 
					score_db=self.scores,
					mutationOpts={
						'degree': self._mutationDegree,
						'library': self._library,
					}, 
					refinementOpts=self._refinementOpts, 
					optimizationOpts=self._optimizationOpts, 
					retain_models=self._retain_models
				) for i in range(self._nChildren)
			)
		return

	def _findSurvivors(self, generation_index, T=400.):
		if generation_index < 0:
			return self.parent
		else:
			scores = self.scores[generation_index, :]
			p = boltzmann_p(scores, T=400.)
			order = np.argsort(p)[::-1]
			total = np.asarray([np.sum(p[0:i+1]) for i in order], dtype=int)
			last_order_idx = np.argwhere(total >= self._survivalCutoff)[0][0]
			survivors = [order[i] for i in range(last_order_idx+1)]
			survivors = [mutant_filename_constructor(self._workdir, generation_index, x) for x in survivors]
			return survivors

	def run(self, verbose=0):
		# if reset:
		# 	self.survivors = [self.parent]
		for i in range(self._nGenerations):
			if verbose != 0:
				print('*** Protean:Evolve - Generating and Evaluating Children for Generation %d ***' % i+1)
			self._generateChildren(generation_index=i)
			# if verbose:
			# 	print('*** Protean:Evolve - Identifying Survivors for Generation %d ***' % i)
			# survivors = self._findSurvivors(generation_index=i)
		if verbose != 0:
			print('*** Protean:Evolve - Genetic Algorithm Complete! ***')
		self._notRun = False
		return

