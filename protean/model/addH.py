from protean.optimize.defaults import _ForceFieldKernel_ImplicitSolvent
from protean.optimize.EnergyMinimization import EnergyMinimization

from simtk.openmm import app


def protonate_protein(topology, positions, variants=None, forcefield=None, platform=None, pH=7.0,
	minimize=False, **kwargs):
	"""
	Summary:
	========
	Use OpenMM modeller method to add back in hydrogens to structures

	See Modeller.addHydrogens for description of variants, forcefield, and platform
	"""
	if forcefield is None:
		forcefield = _ForceFieldKernel_ImplicitSolvent()

	modeller = app.Modeller(topology, positions)
	protonationStates = modeller.addHydrogens(forcefield, pH=pH, platform=platform, variants=variants)

	newtop = modeller.topology
	newpos = modeller.positions

	if minimize:
		optimize = EnergyMinimization(newtop, newpos, **kwargs)
		optimize.run()
		newtop = optimize.getTopology()
		newpos = optimize.getPositions()

	return newtop, newpos