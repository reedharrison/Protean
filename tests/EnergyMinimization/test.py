import os as os
import sys as sys
sys.path.append('../../')

from protean.optimize.EnergyMinimization import EnergyMinimization
from protean.utilities.conversion import openmm2mdtraj
from simtk.openmm import app, Platform

pdbfile = 'barnase_barstar_addH.pdb'

pdb = app.PDBFile(pdbfile)

platform = Platform.getPlatformByName('CUDA')

optimizer = EnergyMinimization(pdb.topology, pdb.positions, platform=platform,
	maxIterations=10000)

optimizer.run()

newtrj = openmm2mdtraj(optimizer.getTopology(), optimizer.getPositions())
newtrj.save('minimized.pdb')
