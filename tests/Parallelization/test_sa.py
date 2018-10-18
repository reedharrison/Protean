from protean.utilities.conversion import openmm2mdtraj, mdtraj2openmm
from protean.optimize.SimulatedAnnealing import SimulatedAnnealing

def optimize(trj, moleculeSelections=['chainid 0 1', 'chainid 0', 'chainid 1']):
	topology, positions = mdtraj2openmm(trj)
	optimizer = SimulatedAnnealing(topology=topology, positions=positions, moleculeSelections=moleculeSelections)
	optimizer.run()
	newtrj = openmm2mdtraj(optimizer.getTopology(), optimizer.getPositions())
	score = optimizer.getScore()
	return (newtrj, score)
