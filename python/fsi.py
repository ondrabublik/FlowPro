import sys
import os

import flowpro as fp

structureSolverPath = os.path.join(fp.flowProPath, 'toolbox', 'elasticity')
sys.path.insert(1, structureSolverPath)


from structureServer import structureServer
from structureSolver import StructureSolver


def run():
	meshDirName = 'turekhron'

	rho = 2
	nu = 0.4
	E = 1e6
	# E = 5.6e6

	pRef = 1e5
	rhoRef = 1.23
	lRef = 1

	fp.run()

	currentPath = os.getcwd()
	os.chdir(structureSolverPath)
	try:
		fp.args('turekhron/dynamics')
		geoPath, simPath, outPath = fp.getPath()
		# meshDirPath = os.path.join(structureSolverPath, 'mesh', meshDirName)
		solver = StructureSolver(geoPath, simPath, outPath)
		structureServer('localhost', 5767, solver, pRef, rhoRef, lRef)
	finally:
		os.chdir(currentPath)


if __name__ == '__main__':
	run()
