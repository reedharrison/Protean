import numpy as np
import simtk.unit as unit

AVO = unit.AVOGADRO_CONSTANT_NA.value_in_unit(unit.mole**-1)
KB = unit.constants.BOLTZMANN_CONSTANT_kB.value_in_unit(unit.kilojoule / unit.kelvin)


def boltzmann_f(x, T=300.):
	if isinstance(x, list):
		x = np.array(x, dtype=float)
	return np.exp(-1.*x / (AVO*KB*T))

def boltzmann_p(scores, T=400.):
	# Note that if T (temperature) is too low, then this form may be unstable
	# temperature may need to be increased to non-physiological values to 
	# prevent Nan's. There may be a more numerically stable solution to this
	# problem.
	# MAX_FLOAT64 = 1.7976931348623157e+308 # cannot represent a number larger than this, force values larger to this particular value.
	if isinstance(scores, list):
		scores = np.array(scores, dtype=float)
	q = boltzmann_f(scores, T=T)
	if np.isinf(q.sum()):
		raise Exception('Overflow in protean.evolve.fitness.boltzmann_p! Try setting T to a larger value.')
	p = q / q.sum()
	return p
