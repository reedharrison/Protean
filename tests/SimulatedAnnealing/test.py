import os as os
import sys as sys
sys.path.append('../../')

from protean.optimize.SimulatedAnnealing import SimulatedAnnealing
from protean.optimize.defaults import _IntegratorKernel_LangevinIntegrator, _ForceFieldKernel_ImplicitSolvent
from protean.utilities.conversion import openmm2mdtraj
from simtk.openmm import app, Platform

pdbfile = 'barnase_barstar_amber.pdb'
prefix = 'barnase_barstar'

pdb = app.PDBFile(pdbfile)
ff = _ForceFieldKernel_ImplicitSolvent()
modeller = app.Modeller(pdb.topology, pdb.positions)
protonationStates = modeller.addHydrogens(ff, pH=7.4)

platform = Platform.getPlatformByName('CUDA')

# def Integrator(temperature):
# 	"""
# 	|      Parameters
# 	|      ----------
# 	|      K_r : integer, default: 2
# 	|          Number of geodesic drift steps.
# 	|
# 	|      temperature : np.unit.Quantity compatible with kelvin, default: 298.0*unit.kelvin
# 	|         Fictitious "bath" temperature
# 	|
# 	|      collision_rate : np.unit.Quantity compatible with 1/picoseconds, default: 1.0/unit.picoseconds
# 	|         Collision rate
# 	|
# 	|      timestep : np.unit.Quantity compatible with femtoseconds, default: 1.0*unit.femtoseconds
# 	|         Integration timestep
# 	|
# 	|      constraint_tolerance : float, default: 1.0e-8
# 	|          Tolerance for constraint solver
# 	|
# 	|      measure_shadow_work : boolean, default: False
# 	|          Accumulate the shadow work performed by the symplectic substeps, in the global `shadow_work`
# 	|
# 	|      measure_heat : boolean, default: False
# 	|          Accumulate the heat exchanged with the bath in each step, in the global `heat`
# 	|
# 	|      References
# 	|      ----------
# 	|      [Leimkuhler and Matthews, 2016] Efficient molecular dynamics using geodesic integration and solvent-solute splitting
# 	|      http://rspa.royalsocietypublishing.org/content/472/2189/20160138
# 	"""
# 	from simtk import unit
# 	# import simtk.openmm as mm
# 	from openmmtools.integrators import GeodesicBAOABIntegrator

# 	K_r = 2
# 	collision_rate = 1/unit.picoseconds
# 	constraint_tolerance = 1.0e-8
# 	timestep = 2.0*unit.femtoseconds

# 	integrator = GeodesicBAOABIntegrator(
# 		K_r=K_r, 
# 		temperature=temperature,
# 		constraint_tolerance=constraint_tolerance,
# 		timestep=timestep
# 	)
# 	return integrator

optimizer = SimulatedAnnealing(topology=modeller.topology, positions=modeller.positions, platform=platform, 
	# constSimulationTemp=310.,
	# restrainAtoms='protein and name CA',
	moleculeSelections=['chainid 0 1', 'chainid 0', 'chainid 1'],
	# temperatureRange=(330., 270.),
	# simulationSteps=500,
	# temperatureSteps=(3, 10),
	# constructIntegrator=Integrator
	)

optimizer.minimize()

newtrj = openmm2mdtraj(optimizer.getTopology(), optimizer.getPositions())
newtrj.save(prefix + '_addH.pdb')

optimizer.run(verbose=True)

newtrj = openmm2mdtraj(optimizer.getTopology(), optimizer.getPositions())
newtrj.save(prefix + '_SAout.pdb')
