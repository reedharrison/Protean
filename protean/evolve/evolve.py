from __future__ import print_function, division

import os as os
import numpy as np
import pickle as p
# import zarr as zarr
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
	def __init__(self, parent, workdir=None, nGenerations=3, nChildren=100, percentRecombination=0.05,
		mutationDegree=1, propogateParent=True, survivalCutoff=0.2, nThreads=-2, library=None,
		sites=None, selection='protein', atom_indices=None, retain_models=True,
		refinementOpts=None, optimizationOpts=None, boltzmannFactor=10.):
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
		self._percentRecombination = percentRecombination
		self._mutationDegree = int(mutationDegree)
		self._propogateParent = bool(propogateParent)
		self._survivalCutoff = float(survivalCutoff)
		self._library = library
		self._nThreads = nThreads
		self._retain_models = retain_models
		self._boltzmannFactor = boltzmannFactor

		if atom_indices is not None:
			self._sites = generate_sites(pdb, indices=atom_indices)
		else:
			self._sites = generate_sites(pdb, selection=selection)

		if refinementOpts is None:
			self._refinementOpts = {'variants':None, 'forcefield':None, 'platform':None, 'pH':7.0, 'minimize':True}
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
		self._checkpoint = os.path.join(self._workdir, 'evolve.cpt') # checkpoint file for to dump results with pickle

		# Persistence mode: ‘r’ means read only (must exist); ‘r+’ means read/write (must exist); 
		#	‘a’ means read/write (create if doesn’t exist); ‘w’ means create (overwrite if exists); 
		#	‘w-‘ means create (fail if exists).

		# n_chunks = np.min([self._nChildren, abs(self._nThreads) * 10])

		self.scores = np.zeros((self._nGenerations, self._nChildren), dtype=float)
		# scoreSynchronizer = zarr.ProcessSynchronizer(os.path.join(self._workdir, 'scores.sync'))
		# self.scores = zarr.open(os.path.join(self._workdir, 'scores.zarr'),
		# 	mode='w', 
		# 	shape=(self._nGenerations, self._nChildren),
		# 	chunks=(1, self._nChildren),
		# 	# chunks=(1, n_chunks),
		# 	dtype=float,
		# 	# synchronizer=scoreSynchronizer
		# 	)

		self.sequences = np.empty((self._nGenerations, self._nChildren), dtype=object)
		# seqSynchronizer = zarr.ProcessSynchronizer(os.path.join(self._workdir, 'sequences.sync'))
		# self.sequences = zarr.open(os.path.join(self._workdir, 'sequences.zarr'),
		# 	mode='w',
		# 	shape=(self._nGenerations, self._nChildren),
		# 	chunks=(1, self._nChildren),
		# 	# chunks=(1, n_chunks),
		# 	dtype=str,
		# 	# synchronizer=seqSynchronizer
		# 	)

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

	# def _generateStructures(self, indices, templates=None, verbose=0, **kwargs):
	# 	# hidden method to generate structure from sequence.
	# 	# useful if structural file gets deleted. saves to work
	# 	# directory for persistence.
	# 	#
	# 	# indicies is a tuple of the form (gen idx, child idx)
	# 	sequences = [self.sequences[g, c] for g, c in indices]
	# 	if templates is None:
	# 		templates = self.parent

	# 	trj = None
	# 	with Parallel(n_jobs=self._nThreads, prefer='processes', verbose=verbose) as parallel:
	# 		trj = parallel(delayed(_generateStructure)(sequence=x, templates=templates, **kwargs) for x in sequences)
	# 	return trj

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

	def _generateChildren(self, generation_index, verbose=0, firstChild=0):
		# parallel generation of children within a single generation
		# call recombine survivors here since natural variation should not
		# replace structures from recombination...
		firstParent = self._propogateParent
		if firstParent: firstParent = self.parent

		with Parallel(backend='loky', n_jobs=self._nThreads, verbose=verbose) as parallel:
			results = parallel(delayed(_variationKernel)
				(
					filenames=self._findSurvivors(generation_index-1),
					workdir=self._workdir,
					generation_index=generation_index, 
					child_index=i, 
					sites=self._sites, 
					# seq_db=self.sequences, 
					# score_db=self.scores,
					mutationOpts={
						'degree': self._mutationDegree,
						'library': self._library,
					}, 
					refinementOpts=self._refinementOpts, 
					optimizationOpts=self._optimizationOpts, 
					retain_models=self._retain_models,
					percentRecombination=self._percentRecombination,
					firstParent=firstParent
				) for i in range(firstChild, self._nChildren)
			)

			scores = []
			sequences = []
			for result in results:
				scores.append(result[0])
				sequences.append(result[1])
			self.scores[generation_index, firstChild:] = scores
			self.sequences[generation_index, firstChild:] = sequences

		with open(self._checkpoint, 'wb') as h:
			p.dump(self, h)

		return

	def _findSurvivors(self, generation_index):
		k = self._boltzmannFactor
		if generation_index < 0:
			return [self.parent]
		elif k > 0.:
			try:
				scores = self.scores[generation_index, :]
				p = boltzmann_p(scores, k=k)
				order = np.argsort(p)[::-1]
				total = np.asarray([np.sum(p[0:i+1]) for i in order], dtype=float)
				last_order_idx = np.argwhere(total >= self._survivalCutoff)[0][0]
				survivors = [order[i] for i in range(last_order_idx+1)]
				survivors = [mutant_filename_constructor(self._workdir, generation_index, x) for x in survivors]
				return survivors
			except:
				print('WARNING: could not calculate Boltzmann probabilities, only propagating most fit sequence')
				scores = self.scores[generation_index, :]
				indmin = np.argmin(scores)
				return [mutant_filename_constructor(self._workdir, generation_index, indmin)]
		elif k <= 0.:
			scores = self.scores[generation_index, :]
			indmin = np.argmin(scores)
			return [mutant_filename_constructor(self._workdir, generation_index, indmin)]

	def run(self, verbose=0):
		# if reset:
		# 	self.survivors = [self.parent]
		for i in range(self._nGenerations):
			if verbose != 0:
				print('*** Protean:Evolve - Generating and Evaluating Children for Generation %d ***' % (i))
			self._generateChildren(generation_index=i)
			# if verbose:
			# 	print('*** Protean:Evolve - Identifying Survivors for Generation %d ***' % i)
			# survivors = self._findSurvivors(generation_index=i)
		if verbose != 0:
			print('*** Protean:Evolve - Genetic Algorithm Complete! ***')
		self._notRun = False
		return

	def restart(self, verbose=0, generation_index=None, child_index=None):
		if generation_index is None:
			genMask = [any([x is None for x in self.sequences[i,:]]) for i in range(self._nGenerations)]
			genIdx = [i for i, flag in enumerate(genMask) if flag][0]
		else:
			genIdx = generation_index

		if child_index is None:
			childMask = [self.sequences[genIdx, j] is None for j in range(self._nChildren)]
			childIdx = [i for i, flag in enumerate(childMask) if flag][0]
		else:
			childIdx = child_index

		for i in range(genIdx, self._nGenerations):
			if verbose != 0:
				print('*** Protean:Evolve - Generating and Evaluating Children for Generation %d ***' % (i))
			if i == genIdx:
				self._generateChildren(generation_index=i, firstChild=childIdx)
			else:
				self._generateChildren(generation_index=i)
		if verbose != 0:
			print('*** Protean:Evolve - Genetic Algorithm Complete! ***')
		self._notRun = False
		return

	def p(self):
		p = boltzmann_p(self.scores[:,:], k=self._boltzmannFactor)
		return p

	def rank(self, n=0):
		# p = self.p(T=T)
		scores = self.scores
		order = np.argsort(scores, axis=None)#[::-1]
		indices = [(i, j) for i, j in zip(*np.unravel_index(order, dims=scores.shape))]
		if n <= 0:
			return indices
		elif n > 0:
			return indices[0:n]
