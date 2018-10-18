"""
Generate Homology Model
=======================
"""

import os as os
import sys as sys
sys.path.append('../../')

from protean.model.homology import homology_model
from protean.model.addH import protonate_protein
from protean.utilities.conversion import mdtraj2openmm, openmm2mdtraj

import mdtraj as md
from simtk.openmm import app


sequence = 'VSDVPRDLEVVAATPTSLLISWDAPAVTVRYYRITYGETGGNSPVQEFTVPGSKSTATISGLKPGVDYTITVYADSFEDPTCTLVTSSGAFQYWISINYRT'
template = md.load('model_PDL1-Mb_final.pdb')
template = template.atom_slice(template.top.select('chainid 0'))

model = homology_model(sequence, templates=[template])
topology, positions = mdtraj2openmm(model)

topology, positions = protonate_protein(topology, positions, minimize=True, maxIterations=1000)
openmm2mdtraj(topology, positions).save('Mb.pdb')


"""
Add back in PDL1 from chain B
=============================
"""

oldtrj = md.load('model_PDL1-Mb_final_addOXT.pdb')
newtrj = md.load('Mb.pdb')

chainA_ref_idx = oldtrj.top.select('chainid 0 and not element H')
chainB_idx = oldtrj.top.select('chainid 1 and not element H')

chainA_ref = oldtrj.atom_slice(chainA_ref_idx)
chainB = oldtrj.atom_slice(chainB_idx)

chainA_idx = newtrj.top.select('chainid 0 and name CA')
chainA_ref_idx = oldtrj.top.select('chainid 0 and name CA')
chainB_idx = oldtrj.top.select('chainid 1 and name CA')

chainA = newtrj.superpose(oldtrj, atom_indices=chainA_idx , ref_atom_indices=chainA_ref_idx)
chainA_idx = newtrj.top.select('chainid 0 and not element H')
chainA = chainA.atom_slice(chainA_idx)

trj = chainA.stack(chainB)
trj.save('Mb-PDL1.pdb')


"""
Minimize and perform short MD
=============================
"""

from protean.optimize.defaults import calcKappa
from protean.optimize.defaults import _ForceFieldKernel_ImplicitSolvent
from protean.utilities.conversion import openmm2mdtraj

from simtk import unit
from simtk.openmm import app, Platform
import simtk.openmm as mm

from openmmtools.integrators import GeodesicBAOABIntegrator

def getPositions(context):
	positions = context.getState(getPositions=True).getPositions()
	return positions

pdbfile = 'Mb-PDL1.pdb'
prefix = 'PDL1-Mb'

temperature = 300 * unit.kelvin
timestep = int(2)
n_steps = int(100000) * int(4)

platform = Platform.getPlatformByName('CUDA')
forcefield = _ForceFieldKernel_ImplicitSolvent()
integrator = GeodesicBAOABIntegrator(K_r=2, timestep=timestep*unit.femtoseconds, temperature=temperature)

pdb = app.PDBFile(pdbfile)
modeller = app.Modeller(pdb.topology, pdb.positions)
protonationStates = modeller.addHydrogens(forcefield, pH=7.0)

system = forcefield.createSystem(modeller.topology,
	nonbondedMethod=app.CutoffNonPeriodic,
	nonbondedCutoff=1.2*unit.nanometers,
	implicitSolvent='GBn2',
	soluteDielectric=20,
	solventDielectric=78.5,
	implicitSolventKappa=calcKappa(
		temperature=temperature.value_in_unit(unit.kelvin),
		dielectric=78.5,
		ion=0.15,
		z=1.
	),
	constraints=app.HBonds
)

simulation = app.Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

simulation.minimizeEnergy(maxIterations=100000, tolerance=2.5*unit.kilojoule_per_mole)

openmm2mdtraj(topology=modeller.topology, positions=getPositions(simulation.context)).save(prefix + '_minimized.pdb')

simulation.context.setVelocitiesToTemperature(temperature)

simulation.reporters.append(app.DCDReporter(prefix + '_trajectory.dcd', 100*timestep))
simulation.reporters.append(app.StateDataReporter(sys.stdout, 100*timestep, step=True, 
	potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
	speed=True, totalSteps=n_steps, separator='\t'))

simulation.step(n_steps)



# optimizer = SimulatedAnnealing(topology=modeller.topology, positions=modeller.positions,
# 	moleculeSelections=['chainid 0 1', 'chainid 0', 'chainid 1'])

# optimizer.minimize(tolerance=2.5*unit.kilojoule_per_mole, maxIterations=1000000)
# newtrj = openmm2mdtraj(optimizer.getTopology(), optimizer.getPositions())
# newtrj.save(prefix + '_addH.pdb')

# optimizer.run(verbose=True)
# newtrj = openmm2mdtraj(optimizer.getTopology(), optimizer.getPositions())
# newtrj.save(prefix + '_SAout.pdb')


# sequence = 'AVSDVPRKLEVVAATPTSLLISWDAPTMVTQYYRITYGETGGNSPVQEFTVPGSKSTATISGLKPGVDYTITVYAVSFEDPTCTLVTSSGAFQYWISINYRTS'

# ref1 = md.load('1ttg.pdb') # Mb
# ref2 = md.load('5dc4.pdb') # Mb chain B
# ref3 = md.load('5jds.pdb') # Nb chain B

# atms2 = ref2.top.select('chainid 1 and protein')
# ref2 = ref2.atom_slice(atms2)
# ref2.save('5dc4_chainB.pdb')

# atms3 = ref3.top.select('chainid 1 and protein')
# ref3 = ref3.atom_slice(atms3)
# ref3.save('5jds_chainB.pdb')

# templates = [ref1, ref2]

# model = homology_model(sequence, templates=templates)
# topology, positions = mdtraj2openmm(model)

# topology, positions = protonate_protein(topology, positions, minimize=True, maxIterations=1000)
# openmm2mdtraj(topology, positions).save('PDL1-Mb_NbInsert.pdb')
