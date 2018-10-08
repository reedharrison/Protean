from __future__ import print_function, division

import os as os
import tempfile as tf
from random import randint

import numpy as np
import mdtraj as md

from modeller import environ, model, alignment, selection
from modeller.automodel import automodel, assess
from modeller.scripts import complete_pdb

from protean.utilities.conversion import model_to_trj, trj_to_models


def create_env(seed=None, hydrogen=False, 
	atom_files_directory=['../atom_files'],
	topology='$(LIB)/top_heav.lib',
	parameters='$(LIB)/par.lib'):
	"""
	Summary:
	========
		Create a modeller environment. This is
		a convenience function that represents 
		our preferred internal configuration
		for homology modelling.
	"""
	if seed is None:
		seed = randint(-50000, 2)

	env = environ(rand_seed=seed) # between 2 and -50000 to not use default rng
	env.io.hydrogen = hydrogen
	env.io.atom_files_directory = atom_files_directory
	env.libs.topology.read(file=topology)
	env.libs.parameters.read(file=parameters)

	return env

def homology_model(sequence, templates, env=None, 
	n_models=1, return_all_models=False):

	prevdir = os.getcwd()
	
	def load_pdbs(pdblist):
		pdb = md.load(pdblist[0])
		for pdbfile in pdblist[1:]:
			pdb = pdb.join(md.load(pdbfile))
		return pdb

	models = []
	template_files = []
	
	with tf.TemporaryDirectory() as tmpdir:
		try:
			# Change to temporary work directory
			os.chdir(os.path.expanduser(tmpdir))

			# Create Modeller environment
			if env is None:
				env = create_env()
			
			# Prepare templates and construct alignment
			aln = alignment(env)
			for i, template in enumerate(templates):
				mdl = None
				if isinstance(template, model):
					# template is Model from Modeller
					mdl = template
				elif isinstance(template, md.Trajectory):
					# template is Trajectory from MDTraj
					mdl = trj_to_models(template, env=env)[-1]
				elif isinstance(template, str):
					# template is filename
					mdl = complete_pdb(filename=template, env=env)
				template_files.append('template%d.pdb' % i)
				mdl.write(template_files[-1])
				aln.append_model(mdl, atom_files=template_files[-1], 
					align_codes=template_files[-1])

			aln.append_sequence(sequence)
			aln[-1].code = 'unknown'
			aln.malign()
			aln.write(file='alignment.pir', alignment_format='PIR')

			# Generate homology models
			a = automodel(env, alnfile='alignment.pir', knowns=template_files, 
				sequence='unknown', assess_methods=(assess.DOPE,))
			a.starting_model = 1
			a.ending_model = n_models
			a.make()

			# (Optional) Retain only best scoring model
			outputs = [x['name'] for x in a.outputs]
			models = load_pdbs(outputs)
			if return_all_models is False:		
				scores = [x['DOPE score'] for x in a.outputs]
				best = np.argmin(scores)
				models = models[best]

			# Return to original working directory
			os.chdir(prevdir)

		# If method fails, return to original working directory
		except:
			print('Error: failed to generate homology model')
			os.chdir(prevdir)

	return models

def build_sequence(sequence, env=None):
	prevdir = os.getcwd()
	with tf.TemporaryDirectory() as tmpdir:
		try:
			# Change to temporary work directory
			os.chdir(os.path.expanduser(tmpdir))
			if env is None:
				env = create_env()
			mdl = model(env)
			mdl.build_sequence(sequence)
			trj = model_to_trj(mdl)
			# Return to original working directory
			os.chdir(prevdir)
		except:
			print('Error: failed to generate homology model')
			os.chdir(prevdir)
	return trj
