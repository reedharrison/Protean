import os as os
import sys as sys

import mdtraj as md
import numpy as np
import pickle as p

import seaborn as sns
import matplotlib.pyplot as plt
plt.ion()

from Bio import SeqIO, Seq, Alphabet, motifs


"""
Methods
=======
"""

def loadEvolver(filename):
	with open(filename, 'rb') as h:
		evolver = p.load(h)
	return evolver

def read_msa(filename, fmt='fasta'):
	msa = []
	with open(filename, 'r') as h:
		for record in SeqIO.parse(h, fmt):
			msa.append(record)
	return msa

def write_msa(records, filename, fmt='fasta'):
	with open(filename, 'w') as h:
		SeqIO.write(records, h, fmt)
	return

def extractEvolvedSequences(evolver, alphabet=Alphabet.IUPAC.IUPACProtein(),
	generation_indices=None, children_indices=None):
	msa = []

	sequences = evolver.sequences
	if generation_indices is not None:
		sequences = sequences[generation_indices,:]
	if children_indices is not None:
		sequences = sequences[:, children_indices]
	sequences = sequences.flatten().tolist()

	for i, sequence in enumerate(sequences):
		gid, cid = np.unravel_index(i, (evolver._nGenerations, evolver._nChildren))
		name = 'mutant_g%d_c%d' % (gid, cid)
		desc = 'child %d from generation %d' % (cid, gid)
		seq = Seq.Seq(sequence, alphabet=alphabet)
		msa.append(SeqIO.SeqRecord(
			seq=seq,
			id=str(i),
			name=name,
			description=desc
			)
		)
	return msa

def reweightedMotif(sequences, fitness, k=10.):
	alphabet = sequences[0].seq.alphabet
	n_records = len(sequences)
	n_positions = len(sequences[0])
	n_letters = len(alphabet.letters)
	assert all([len(x) == n_positions for x in sequences])
	assert all([x.seq.alphabet == alphabet for x in sequences])

	factors = fitness / fitness.min()
	probs = np.exp(factors / float(k))
	probs = probs / probs.sum()
	probs = probs.reshape((-1, 1, 1))
	probs = np.repeat(probs, n_positions, axis=1)
	probs = np.repeat(probs, n_letters, axis=2)

	counts = np.zeros((n_records, n_positions, n_letters))

	for i, record in enumerate(sequences):
		for j in range(n_positions):
			for k in range(n_letters):
				if record[j] == alphabet.letters[k]:
					counts[i, j, k] += 1

	counts = counts * probs
	counts = counts * n_records / counts.sum()
	counts = counts.sum(axis=0)
	d = {}
	for i, letter in enumerate(alphabet.letters):
		d[letter] = counts[:, i].tolist()

	m = motifs.Motif(counts=d, alphabet=alphabet)
	return m

def sequenceLogo(evolver, filename=None, k=10., chainid=0, generation_indices=None, 
	children_indices=None, alphabet=Alphabet.IUPAC.IUPACProtein(), format='SVG'):
	firstPosition = np.min([x[1] for x in evolver._sites if x[0] == chainid])
	lastPosition = np.max([x[1] for x in evolver._sites if x[0] == chainid]) + 1
	sequences = extractEvolvedSequences(evolver, generation_indices=generation_indices,
		children_indices=children_indices, alphabet=alphabet)
	sequences = [x[firstPosition:lastPosition] for x in sequences]

	scores = np.asarray([x for x in evolver.scores.flatten() if not np.isnan(x)])
	sequences = [x for x, y in zip(sequences, evolver.scores.flatten()) if not np.isnan(y)]

	m = reweightedMotif(sequences, fitness=scores, k=10.)
	if filename is None:
		return m
	else:
		m.weblogo(filename, format='SVG', first_index=firstPosition,
			logo_start=firstPosition, logo_end=lastPosition-1)
		return

def plotScoreDistribution(evolver, byGeneration=True, ax=None, **kwargs):
	if ax is None:
		fig, ax = plt.subplots(1,1)
	if byGeneration is True:
		scores = evolver.scores.flatten()
		scores_size = evolver.scores.size
		scores_shape = evolver.scores.shape
		generations = [np.unravel_index(i, scores_shape)[0] for i in range(scores_size)]
		h = sns.violinplot(x=generations, y=scores, ax=ax, **kwargs)
		ax.set_xlabel('Generation index')
		ax.set_ylabel('Score (kJ/mol)')
	elif byGeneration is False:
		scores = evolver.scores.flatten()
		h = sns.distplot(scores, ax=ax, **kwargs)
		ax.set_xlabel('Score (kJ/mol)')
		ax.set_ylabel('Density')
	if ax is None:
		fig.tight_layout()
		return fig, ax
	else:
		return h
