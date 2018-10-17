from protean.optimize.defaults import _ForceFieldKernel_ImplicitSolvent
from protean.optimize.defaults import _SystemKernel_ImplicitSolvent
from protean.optimize.defaults import _IntegratorKernel_SimulatedAnnealing
from protean.optimize.defaults import _boltzmannProbability

from protean.utilities.conversion import openmm2mdtraj

from random import random

from simtk import unit as u
from simtk.openmm import app
import simtk.openmm as mm

from openmmtools.integrators import GeodesicBAOABIntegrator

import numpy as np


def _getPositions(context):
	positions = context.getState(getPositions=True).getPositions()
	return positions

class SimulatedAnnealing:
	"""
	Note: immobile atoms doesn't work right now, and I'm not sure why. Unconstrained atoms react unpredictably. Restrained atoms are fine though
	"""
	def __init__(self, topology, positions, temperatureRange=(350., 250.), temperatureSteps=(11, 50), 
		simulationSteps=100, forcefield=None, constructIntegrator=None, constructSystem=None, platform=None,
		moleculeSelections=['protein'], constSimulationTemp=None, constrainAtoms=None, restrainAtoms=None,
		restraintConstant=5.0*u.kilojoules/(u.angstrom**2)):

		self._temperatureRange = temperatureRange # (Ti, Tf)
		self._temperatureSteps = temperatureSteps # (number of temps, number of runs at each temp)
		self._simulationSteps = simulationSteps   # number of steps for integrator
		self._constSimulationTemp = constSimulationTemp

		initialTemp = constSimulationTemp # Simulation temperature is constant if this is not None in which case only accept/reject is affected by temperature)
		if initialTemp is None:
			initialTemp = np.max(temperatureRange)

		if forcefield is None:
			forcefield = _ForceFieldKernel_ImplicitSolvent()
		self._forcefield = forcefield

		self._platform = platform

		if constructIntegrator is None:
			constructIntegrator = _IntegratorKernel_SimulatedAnnealing
		self._constructIntegrator = constructIntegrator

		if constructSystem is None:
			constructSystem = _SystemKernel_ImplicitSolvent
		self._constructSystem = constructSystem

		# System Containers
		self._keys = []
		self._atomMap = {}
		self._topology = {}
		self._positions = {}
		self._simulation = {}
		self._energy = {}

		# Partition System
		trj = openmm2mdtraj(topology=topology, positions=positions)
		for i, selection in enumerate(moleculeSelections):
			self._keys.append(i)
			self._atomMap[i] = trj.top.select(selection)
			subset = trj.atom_slice(self._atomMap[i])
			self._topology[i] = subset.top.to_openmm()
			self._positions[i] = subset.openmm_positions(-1)

		# Configure Constraints and Restraints
		if isinstance(constrainAtoms, str):
			self._constrainAtoms = trj.top.select(constrainAtoms)
		else:
			self._constrainAtoms = constrainAtoms
		if isinstance(restrainAtoms, str):
			self._restrainAtoms = trj.top.select(restrainAtoms)
		else:
			self._restrainAtoms = restrainAtoms
		self._restraintConstant = restraintConstant

		# Initialize simulation objects
		self._initializeSimulations(temperature=initialTemp)
		# self._initializeEnergies(temperature=initialTemp)
		for i in self._keys:
			self._energy[i] = self._calculateEnergy(self._simulation[i])
		return

	def __str__(self):
		dE = self.getScore()
		line = 'Simulated Annealing Optimization Score = ' + dE.format('%0.2f')
		return line
	
	def _initializeSimulations(self, temperature):
		for i in self._keys:
			integrator = self._constructIntegrator(temperature=temperature)
			system = self._constructSystem(topology=self._topology[i],
				forcefield=self._forcefield, temperature=temperature)
			if self._constrainAtoms is not None:
				self._applyConstraints(system, key=i)
			if self._restrainAtoms is not None:
				self._applyRestraints(system, key=i)
			self._simulation[i] = app.Simulation(self._topology[i], system=system, 
				integrator=integrator, platform=self._platform)
			self._simulation[i].context.setPositions(self._positions[i])
		return

	# def _initializeEnergies(self, temperature):
	# 	for i in self._keys:
	# 		self._simulation[i].context.setVelocitiesToTemperature(temperature)
	# 	energies = [self._calculateEnergy(self._simulation[i]) for i in self._keys]
	# 	self._storeEnergy(keys=self._keys, energy=energies)
	# 	return

	def _applyConstraints(self, system, key):
		for idx, atom in zip(self._atomMap[key], self._topology[key].atoms()):
			if idx in self._constrainAtoms:
				system.setParticleMass(atom.index, 0.0)
		return

	def _applyRestraints(self, system, key):
		force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
		force.addGlobalParameter("k", self._restraintConstant)
		force.addPerParticleParameter("x0")
		force.addPerParticleParameter("y0")
		force.addPerParticleParameter("z0")
		for idx, atom, xyz in zip(self._atomMap[key], self._topology[key].atoms(), self._positions[key]):
			if idx in self._restrainAtoms:
				force.addParticle(atom.index, xyz.value_in_unit(u.nanometers))
		system.addForce(force)
		return

	def minimize(self, maxIterations=1000, tolerance=5*u.kilojoule_per_mole):
		key = self._keys[0]
		self._simulation[key].minimizeEnergy(maxIterations=maxIterations, tolerance=tolerance)
		self._updatePositions()
		for i in self._keys:
			self._energy[i] = self._calculateEnergy(self._simulation[i])
		return

	def _calculateEnergy(self, simulation, positions=None):
		if positions is not None:
			simulation.context.setPositions(positions)
		state = simulation.context.getState(getEnergy=True)
		return state.getPotentialEnergy()

	def _subsetPositions(self, positions):
		keys = self._keys
		subset = {i: [] for i in keys}
		for idx, pos in zip(self._atomMap[keys[0]], positions):
			for i in keys[1:]:
				if idx in self._atomMap[i]:
					subset[i].append(pos)
		subset[keys[0]] = positions
		return subset

	def _updatePositions(self):
		key = self._keys[0]
		initialPositions = _getPositions(self._simulation[key].context)
		finalPositions = self._subsetPositions(initialPositions)
		for i in self._keys:
			self._simulation[i].context.setPositions(finalPositions[i])
		return

	def _validateBondLengths(self, cutoffBondLength=0.25*u.nanometers):
		# i = 0 # only validate bonds for largest structure, should be index 0...
		for i in self._keys:
			top = self.getTopology(key=i)
			pos = _getPositions(self._simulation[i].context)
			for a1, a2 in top.bonds():
				x1, y1, z1 = pos[a1.index]
				x2, y2, z2 = pos[a2.index]
				dist = ((x1 - x2)**2. + (y1 - y2)**2. + (z1 - z2)**2.)**0.5
				if dist > cutoffBondLength:
					# print("\n******INVALID BOND LENGTH DETECTED******\n")
					return False
		return True

	def _storePositions(self):
		for i in self._keys:
			self._positions[i] = _getPositions(self._simulation[i].context)
		return

	def _storeEnergy(self, keys, energy):
		for key, val in zip(keys, energy):
			self._energy[key] = val
		return

	def _step(self, temperature):
		nKeys = len(self._keys)
		key = self._keys[0]
		nSteps = self._simulationSteps

		if self._constSimulationTemp is not None:
			self._simulation[key].context.setVelocitiesToTemperature(self._constSimulationTemp)
		else:
			self._simulation[key].context.setVelocitiesToTemperature(temperature)

		self._simulation[key].step(nSteps)
		self._updatePositions()

		newEnergy = [self._calculateEnergy(self._simulation[i]) for i in self._keys]
		oldEnergy = [self._energy[i] for i in self._keys]

		Ei = oldEnergy[0]
		Ef = newEnergy[0]
		if nKeys > 1:
			for val in oldEnergy[1:]:
				Ei -= val
			for val in newEnergy[1:]:
				Ef -= val
		prob = _boltzmannProbability(Ei, Ef, temperature=temperature)

		valid = self._validateBondLengths() # only add structure if bonds are somewhat valid...

		if (random() <= prob) and (valid is True):
			self._storePositions()
			self._storeEnergy(keys=self._keys, energy=newEnergy)
		return

	def run(self, verbose=False, reinitialize=True, maxNoImprovement=5):
		reinitialize = reinitialize & (self._constSimulationTemp is None) # only reinitialize if not constSimulationTemp

		if verbose: print(str(self))

		Ti = np.max(self._temperatureRange)
		Tf = np.min(self._temperatureRange)
		nT, rT = self._temperatureSteps

		counter = 0 # count number iter no improvement
		lastScore = self.getScore()

		for T in np.linspace(Ti, Tf, nT, dtype=float):
			if (T != Ti) and reinitialize:
				self._initializeSimulations(T)
			for i in range(rT):
				self._step(temperature=T)

				if verbose: print(str(self))

				if maxNoImprovement > 0:
					currentScore = self.getScore()

					if  currentScore >= lastScore:
						counter += 1
					else:
						counter = 0

					lastScore = currentScore

					if counter >= maxNoImprovement:
						counter = 0
						break
		return

	def getPositions(self, key=None):
		if key is None:
			key = self._keys[0]
		return self._positions[key]

	def getTopology(self, key=None):
		if key is None:
			key = self._keys[0]
		return self._topology[key]

	def getScore(self):
		energy = [self._energy[i] for i in self._keys]
		dE = energy[0]
		if len(energy) > 1:
			for val in energy[1:]:
				dE -= val
		return dE
