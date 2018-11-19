from protean.optimize.defaults import _ForceFieldKernel_ImplicitSolvent
from protean.optimize.defaults import _SystemKernel_ImplicitSolvent
from protean.optimize.defaults import _IntegratorKernel_EnergyMinimization

from protean.utilities.conversion import openmm2mdtraj

import simtk.unit as unit
from simtk.openmm import app
import simtk.openmm as mm


def _getPositions(context):
	positions = context.getState(getPositions=True).getPositions()
	return positions

class EnergyMinimization:
	def __init__(self, topology, positions, energyTol=5*unit.kilojoule_per_mole, maxIterations=1000,
		forcefield=None, constructIntegrator=None, constructSystem=None, 
		platform=None, restrainAtoms=None, temperature=300.,
		restraintConstant=5.0*unit.kilojoules/(unit.angstrom**2)):

		self._topology = topology
		self._positions = positions
		self._energyTol = energyTol
		self._maxIterations = int(maxIterations)
		self._temperature = temperature

		if forcefield is None:
			forcefield = _ForceFieldKernel_ImplicitSolvent()
		self._forcefield = forcefield

		if constructIntegrator is None:
			constructIntegrator = _IntegratorKernel_EnergyMinimization
		self._constructIntegrator = constructIntegrator

		if constructSystem is None:
			constructSystem = _SystemKernel_ImplicitSolvent
		self._constructSystem = constructSystem

		self._platform = platform

		if restrainAtoms is not None:
			trj = openmm2mdtraj(self._topology, self._positions)
			if isinstance(restrainAtoms, str):
				self._restrainAtoms = trj.top.select(restrainAtoms)
			else:
				self._restrainAtoms = restrainAtoms
		else:
			self._restrainAtoms = restrainAtoms
		self._restraintConstant = restraintConstant

		self._initializeSimulation()
		return

	def _initializeSimulation(self):
		integrator = self._constructIntegrator(temperature=self._temperature*unit.kelvin)
		system = self._constructSystem(topology=self._topology,
			forcefield=self._forcefield, temperature=self._temperature)
		if self._restrainAtoms is not None:
			self._applyRestraints(system)
		self._simulation = app.Simulation(self._topology, system=system, 
			integrator=integrator, platform=self._platform)
		self._simulation.context.setPositions(self._positions)
		return

	def _applyRestraints(self, system):
		force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
		force.addGlobalParameter("k", self._restraintConstant)
		force.addPerParticleParameter("x0")
		force.addPerParticleParameter("y0")
		force.addPerParticleParameter("z0")
		for atom, xyz in zip(self._topology.atoms(), self._positions):
			force.addParticle(atom.index, xyz.value_in_unit(unit.nanometers))
		system.addForce(force)
		return

	def getPositions(self):
		return _getPositions(self._simulation.context)

	def getTopology(self):
		return self._topology

	def run(self):
		self._simulation.minimizeEnergy(maxIterations=self._maxIterations, tolerance=self._energyTol)
		return
		