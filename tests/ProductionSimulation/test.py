import os as os
import sys as sys
sys.path.append('../../')

from protean.optimize.defaults import calcKappa
from protean.optimize.defaults import _ForceFieldKernel_ImplicitSolvent
from protean.utilities.conversion import openmm2mdtraj

from simtk import unit
from simtk.openmm import app, Platform
import simtk.openmm as mm

from openmmtools.integrators import GeodesicBAOABIntegrator

"""
Methods
=======
"""
def getPositions(context):
	positions = context.getState(getPositions=True).getPositions()
	return positions


"""
Setup and Run Simulation
========================
"""
def main():
	pdbfile = 'barnase_barstar_amber.pdb'
	temperature = 300 * unit.kelvin

	timestep = int(4)
	n_steps = int(10000) * int(4)

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

	simulation.minimizeEnergy()

	openmm2mdtraj(topology=modeller.topology, positions=getPositions(simulation.context)).save('minimized.pdb')

	simulation.context.setVelocitiesToTemperature(temperature)

	simulation.reporters.append(app.DCDReporter('trajectory.dcd', 100*timestep))
	simulation.reporters.append(app.StateDataReporter(sys.stdout, 100*timestep, step=True, 
		potentialEnergy=True, temperature=True, progress=True, remainingTime=True, 
		speed=True, totalSteps=n_steps, separator='\t'))

	simulation.step(n_steps)
	return

main()