import os as os
import numpy as np
import mdtraj as md

from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit as u

from math import exp

try:
	from openmmtools.integrators import GeodesicBAOABIntegrator
except:
	print('OpenMMTools not installed, make sure to specify your own integrator for Simulated Annealing')


"""
Globals:
========
"""

# _DEFAULT_FF_SOLUTE = 'amber10.xml'
# _DEFAULT_FF_SOLVENT = 'amber99_obc.xml'
_DEFAULT_FF_SOLUTE = 'amber99sb.xml'
_DEFAULT_FF_SOLVENT = 'amber99_obc.xml'

"""
Common Functions:
=================
"""

def _boltzmannProbability(Ei, Ef, temperature=300.):
	if Ef > Ei:
		T = temperature * u.kelvin
		kb = u.constants.BOLTZMANN_CONSTANT_kB
		avo = u.constants.AVOGADRO_CONSTANT_NA
		return exp((Ei - Ef) / (kb * T * avo))
	else:
		return 1.

def calcKappa(dielectric=78.5, temperature=310., ion=0.15, z=1.):
	from math import sqrt

	vac = 8.854187817620e-12 # Vacuum permittivity [C^2 / J / m]
	kb = u.constants.BOLTZMANN_CONSTANT_kB.value_in_unit(u.joule / u.kelvin)
	avo = u.constants.AVOGADRO_CONSTANT_NA.value_in_unit(u.mole ** -1)
	e = u.constants.elementary_charge.conversion_factor_to(u.coulomb)
	f = u.liter.conversion_factor_to(u.meter ** 3) # convert liter to m^3

	numer = 2.0 * ion * ((z * e) ** 2) * avo 
	denom = dielectric * vac * kb * temperature * f
	kappa = sqrt(numer / denom) / u.meter # kappa is in units of m^-1

	return kappa.in_units_of(u.nanometer ** -1)


"""
Default Kernels:
================
"""

def _ForceFieldKernel_ImplicitSolvent():
	ff = app.ForceField(_DEFAULT_FF_SOLUTE, _DEFAULT_FF_SOLVENT)
	return ff

def _SystemKernel_ImplicitSolvent(topology, forcefield, temperature,
	nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=1.2*u.nanometers,
	implicitSolvent='GBn2', soluteDielectric=20., solventDielectric=78.5,
	ionConcentration=0.15, ionChargeMagnitude=1., constraints='AllBonds'):
	system = forcefield.createSystem(topology,
		nonbondedMethod = nonbondedMethod,
		nonbondedCutoff = nonbondedCutoff,
		implicitSolvent = implicitSolvent,
		soluteDielectric = soluteDielectric,
		solventDielectric = solventDielectric,
		implicitSolventKappa = calcKappa(
			temperature = temperature,
			dielectric = solventDielectric,
			ion = ionConcentration,
			z = ionChargeMagnitude
		),
		constraints = constraints
	)
	return system

# def _IntegratorKernel_SimulatedAnnealing(temperature, 
# 	frictionCoef=1/u.picoseconds, errorTol=0.001, consTol=0.00001):
# 	integrator = mm.VariableLangevinIntegrator(temperature, frictionCoef, errorTol)
# 	integrator.setConstraintTolerance(consTol)
# 	return integrator

def _IntegratorKernel_SimulatedAnnealing(temperature, timestep=2*u.femtoseconds,
	K_r=2, collision_rate=1/u.picoseconds, constraint_tolerance=1e-8):
	integrator = GeodesicBAOABIntegrator(
		temperature = temperature*u.kelvin, 
		K_r = K_r,
		collision_rate = collision_rate,
		timestep = timestep,
		constraint_tolerance = constraint_tolerance
	)
	return integrator

def _IntegratorKernel_EnergyMinimization(errorTol=0.001, **kwargs):
	integrator = mm.VariableVerletIntegrator(errorTol)
	return integrator

def _IntegratorKernel_LangevinIntegrator(temperature, 
	frictionCoef=1/u.picoseconds, timeStep=0.002*u.picoseconds):
	integrator = mm.LangevinIntegrator(temperature, frictionCoef, timeStep)
	return integrator
