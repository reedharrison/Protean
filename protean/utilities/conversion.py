import os as os
import mdtraj as md
import tempfile as tf

from modeller.scripts import complete_pdb

from simtk import unit as u


def mdtraj2openmm(trj, frame=-1):
	return trj.top.to_openmm(), trj.openmm_positions(frame)

def openmm2mdtraj(topology, positions):
	pos = positions.value_in_unit(u.nanometer)
	top = md.Topology().from_openmm(topology)
	return md.Trajectory(pos, topology=top)

def read_model(filename, env=None, transfer_res_num=True):
	if env is None:
		env = create_env()
	model = complete_pdb(env, filename, transfer_res_num=transfer_res_num)
	return model

def model_to_trj(model):
	"""
	Summary:
	========
		Write a model from Modeller to file then parse
		using MDtraj and return trajectory object
	"""
	with tf.TemporaryDirectory() as h:
		tmpdir = os.path.expanduser(h)
		tmpfile = os.path.join(tmpdir, 'mdl.pdb')
		model.write(file=tmpfile)
		trj = md.load(tmpfile)
		# trj.top.create_disulfide_bonds(trj.openmm_positions(-1)) # 10/11/2018, try to update bond block for DISU to prevent topology errors with CYS
	return trj

def trj_to_models(trj, env=None):
	"""
	Summary:
	========
		Convert mdtraj Trajectory to list of Modeller 
		models by saving PDB then reading it back in.

		Note: models will still need patching!
	"""
	if env is None:
		env = create_env()
	models = []
	with tf.TemporaryDirectory() as h:
		tmpdir = os.path.expanduser(h)
		tmpfile = os.path.join(tmpdir, 'mdl.pdb')
		for frame in trj:
			frame = frame.center_coordinates() # added 2018/10/12 to avoid writing errors with out of bound coordinates!
			frame.save(tmpfile)
			model = complete_pdb(env, tmpfile, transfer_res_num=True)
			models.append(model)
	return models
