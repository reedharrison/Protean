from __future__ import print_function, division

import os as os

def mutant_filename_constructor(workdir, generation_index, child_index):
	filename = os.path.join(workdir, 'mutant_g%d_c%d.pdb' % (generation_index, child_index))
	return filename
