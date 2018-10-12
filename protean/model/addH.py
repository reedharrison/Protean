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

	weird_cb_bonds = [(x, y) for x, y in modeller.topology.bonds() if x.name=='CB' and y.name=='CB']
	modeller.delete(weird_cb_bonds) # bonds do not occur between CB atoms ever, remove all such instances (i think modeller does this for CG models or something)

	modeller.topology.createDisulfideBonds(modeller.positions)
	protonationStates = modeller.addHydrogens(forcefield, pH=pH, platform=platform, variants=variants)

	newtop = modeller.topology
	newpos = modeller.positions

	if minimize:
		optimize = EnergyMinimization(newtop, newpos, **kwargs)
		optimize.run()
		newtop = optimize.getTopology()
		newpos = optimize.getPositions()

	return newtop, newpos