from protean.optimize.defaults import _ForceFieldKernel_ImplicitSolvent
from protean.optimize.EnergyMinimization import EnergyMinimization

from simtk.openmm import app
from simtk import unit


def protonate_protein(topology, positions, variants=None, forcefield=None, platform=None, pH=7.0,
	minimize=True, cutoffBondLength=2.5*unit.nanometers, **kwargs):
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

	if not minimize: # minimize if bad bonds detected
		for a1, a2 in newtop.bonds():
			x1, y1, z1 = newpos[a1.index]
			x2, y2, z2 = newpos[a2.index]
			dist = ((x1 - x2)**2. + (y1 - y2)**2. + (z1 - z2)**2.)**0.5
			if dist > cutoffBondLength:
				print('Bond length violations detected:')
				print(a1, a2, dist)
				minimize = True

	if minimize:
		optimize = EnergyMinimization(newtop, newpos, **kwargs)
		optimize.run()
		newtop = optimize.getTopology()
		newpos = optimize.getPositions()

	return newtop, newpos
