from __future__ import print_function, division

import os as os
import numpy as np
import mdtraj as md
import itertools as it

from random import choice
from modeller import model
# from modeller.scripts import complete_pdb

from protean.model.homology import homology_model
from protean.model.addH import protonate_protein
from protean.utilities.conversion import mdtraj2openmm, openmm2mdtraj


"""
Schemes
=======
"""

_DEFAULT_AA_LIBRARY = list('ARNDEQGHILKMFPSTWYV')

def _DEFAULT_MUTATION_SCHEME(type):
	try:
		assert isinstance(type, str)
	except:
		raise Exception('Default mutation scheme only string types')

	aa_set = set('ARNDCEQGHILKMFPSTWYV')
	output = []
	for x in type:
		if x == '-': # this is not strictly necessary to check, but a reminder that non-recognized characters will be stripped...
			output.append('')
		elif x in aa_set:
			output.append(x)
	return ''.join(output)


"""
Methods
=======
"""

def get_sequence(trj, chainSep='/'):
	sequence = []
	if isinstance(trj, model):
		for chain in trj.chains:
			if len(sequence) != 0:
				sequence.append(chainSep)
			for residue in chain.residues:
				sequence.append(residue.code)
	elif isinstance(trj, md.Trajectory):
		for chain in trj.topology.chains:
			if len(sequence) != 0:
				sequence.append(chainSep)
			for residue in chain.residues:
				sequence.append(residue.code)
	elif isinstance(trj, str):
		trj = md.load(trj)
		for chain in trj.chains:
			if len(sequence) != 0:
				sequence.append(chainSep)
			for residue in chain.residues:
				sequence.append(residue.code)
	else:
		raise Exception('No method to get sequence for %s' % type(trj))
	return ''.join(sequence)

def generate_sites(trj, selection=None, indices=None):
	if isinstance(trj, str):
		trj = md.load(trj)
	assert isinstance(trj, md.Trajectory)

	sites = []
	if selection is not None:
		targets = trj.top.select(selection)
		for i, chain in enumerate(trj.top.chains):
			for j, residue in enumerate(chain.residues):
				if any([x.index in targets for x in residue.atoms]): # if any atom from residue is selected then add site!
					sites.append(tuple([i, j]))
	elif indices is not None:
		for i, chain in enumerate(trj.top.chains):
			for j, residue in enumerate(chain.residues):
				if any([x.index in indices for x in residue.atoms]):
					sites.append(tuple([i, j]))
	return sites

# def insert_residues(trj, residues, positions): # I think my current implementation won't require these methods, in growing
# 	"""
# 	"""
# 	return

# def delete_residues():
# 	"""
# 	"""
# 	return

def mutate_sequence(parent, mutations, scheme=None):
	"""
	mutations must be of the form:
		mutations = [{'chainid': x, 'position': y, 'type': z}, ...]
	where each element of the list specifies a single mutation. In the future, deletion
	or deletion could occur at this level by specifying a delete or insert key in type.
	Remember posiion must be an index that starts from 0.

	assume child inherits parent's structure, this may not actually occur...

	if you want to grow sequence, mutations can specify a type that implies expansion such as "RK"
	which would suggest to replace a single position with residues R and K, in that order.

	if you want to shrink a sequence, mutations can not specify '' or some other character not recognized
	by the scheme. This would replace the residue at a desired position with '', removing it.

	caution is warranted when growing/shrinking sequences. it may be better to develop
	some other method for this.
	"""
	newsequence = []
	for i, segment in enumerate(get_sequence(parent).split('/')):
		positions = []
		types = []
		for mutation in mutations:
			if mutation['chainid'] == i:
				positions.append(mutation['position'])
				types.append(mutation['type'])
			order = np.argsort(positions)[::-1] # by adjusting position from right-hand-side, we can avoid perturbing indexing from left...
			positions = [positions[x] for x in order]
			types = [types[x] for x in order]

			for position, type in zip(positions, types):
				segment = apply_mutation_type(segment=segment, position=position, type=type, scheme=scheme)
		newsequence.append(segment)
	return '/'.join(newsequence)

def apply_mutation_type(segment, position, type, scheme=None):
	if scheme is None:
		scheme = _DEFAULT_MUTATION_SCHEME

	residues = list(segment)
	n_positions = len(residues)
	try:
		assert position < n_positions
	except:
		raise Exception('Mutation position does not exist in segment')

	residues[position] = scheme(type)

	return ''.join(residues)

def mutate_structure(parent, mutations, scheme=None, **kwargs):
	child_sequence = mutate_sequence(parent=parent, mutations=mutations, scheme=scheme)
	models = homology_model(sequence=child_sequence, templates=[parent], **kwargs)
	return models

def refine_structures(trjs, **kwargs):
	if isinstance(trjs, list):
		structures = []
		for trj in trjs:
			topology, positions = mdtraj2openmm(trj)
			topology, positions = protonate_protein(topology, positions, **kwargs)
			structures.append(openmm2mdtraj(topology, positions))
		return structures
	else:
		topology, positions = mdtraj2openmm(trjs)
		topology, positions = protonate_protein(topology, positions, **kwargs)
		return openmm2mdtraj(topology, positions)

def random_mutagenesis(parent, sites, degree=1, library=None):
	# sites should be of the form (chain_idx, position_idx)
	if library is None:
		library = _DEFAULT_AA_LIBRARY
	segments = get_sequence(parent, chainSep='/').split('/')
	mutations = []
	for chain_idx, position_idx in choice(list(it.combinations(sites, degree))):
		newtypes = [x for x in library if segments[chain_idx][position_idx] != x]
		mutations.append({
			'chainid': chain_idx, 
			'position': position_idx, 
			'type': choice(newtypes)
		})
	mdl = mutate_structure(parent=parent, mutations=mutations)
	return mdl
