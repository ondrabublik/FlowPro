import sys
import os
import math

import flowpro as fp
import elasticity as ela

structureSolverPath = os.path.join(fp.flowProPath, 'toolbox', 'elasticity')
sys.path.insert(1, structureSolverPath)


from structureServer import structureServer
from structureSolver import StructureSolver


def run():
	# fluidGeom = 'turekhron/dynamics'
	# fluidSim = 'fsi3'
	# fp.args(fluidGeom, fluidSim)
	#
	# structureGeom = 'turekhron/4'
	# structureSim = 'fsi3'
	# ela.args(structureGeom, structureSim)

	geoPath, simPath, outPath = fp.getPath()

	paramDct = fp.getParam()

	model = paramDct['model']
	if model == 'NavierStokesVelocityInlet':
		vRef = float(paramDct['vIn'])
		rhoRef = float(paramDct['rhoIn'])
		pRef = rhoRef * vRef ** 2
	elif model == 'IncompressibleNavierStokes':
		lst = paramDct['VIn'].split(',')
		vRef = 0.
		for str in lst:
			vRef += float(str.strip()) ** 2
		vRef = math.sqrt(vRef)
		rhoRef = float(paramDct['density'])
		pRef = rhoRef * vRef ** 2

	lRef = paramDct.get('lRef')
	if lRef is None:
		lRef = 1
	else:
		lRef = float(lRef)

	structureParamDct = fp.getStructureParam()
	elaGeom = structureParamDct['geometry']

	elaGeomPath = os.path.join(ela.structureSolverPath, 'simulations', elaGeom, 'mesh')

	# mach = float(paramDct['mach'])
	# kappa = float(paramDct['kappa'])
	# pOut = 1 / (mach**2 * kappa)

	fp.run()

	currentPath = os.getcwd()
	os.chdir(structureSolverPath)
	try:
		# fp.args('turekhron/dynamics')
		# geoPath, simPath, outPath = fp.getPath()
		# meshDirPath = os.path.join(structureSolverPath, 'mesh', meshDirName)

		solver = StructureSolver(elaGeomPath, simPath, outPath)
		structureServer('localhost', 5767, solver, pRef, rhoRef, lRef)
	finally:
		os.chdir(currentPath)


if __name__ == '__main__':
	run()
