import sys
import os
import math

import flowpro as fp
import elasticity as ela

structureSolverPath = ela.structureSolverPath
sys.path.insert(1, structureSolverPath)


from structureServer import StructureServer, StructureServerTest, StructureServerALE
from dynamics import Dynamics


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
		kapa = float(paramDct['kappa'])
		mach =float(paramDct['mach'])
		soundSpeed = vRef / mach
		pOut = rhoRef * soundSpeed**2 / kapa
	elif model == 'IncompressibleNavierStokes':
		pOut = float(paramDct['pOut'])
		pOut = 0

	print('outlet pressure = %.3e' % pOut)  # 1.142857142857143e+06

	# model = paramDct['model'] 
	# if model == 'NavierStokesVelocityInlet':
	# 	vRef = float(paramDct['vIn'])
	# 	rhoRef = float(paramDct['rhoIn'])
	# 	pRef = rhoRef * vRef ** 2
	# elif model == 'IncompressibleNavierStokes':
	# 	lst = paramDct['VIn'].split(',')
	# 	vRef = 0.
	# 	for str in lst:
	# 		vRef += float(str.strip()) ** 2
	# 	vRef = math.sqrt(vRef)
	# 	rhoRef = float(paramDct['density'])
	# 	pRef = rhoRef * vRef ** 2

	# lRef = paramDct.get('lRef')
	# if lRef is None:
	# 	lRef = 1
	# else:
	# 	lRef = float(lRef)

	port = paramDct.get('remoteStructureSolverPort')
	if port is None:
		port = 5767
	else:
		port = int(port)

	print('using port %d' % port)

	structureParamDct = fp.getStructureParam()
	elaGeom = structureParamDct['geometry']	

	elaGeomPath = os.path.join(ela.structureSolverPath, 'simulations', elaGeom, 'mesh')

	currentPath = os.getcwd()
	os.chdir(structureSolverPath)
	try:
		animationPath = os.path.join(simPath, 'animation')
		solver = Dynamics(elaGeomPath, simPath, outPath, animationPath)

		regime = structureParamDct.get('regime')
		fp.run()

		print('starting FSI computation')
		server = StructureServer('localhost', port, solver, pOut, regime)
		# if regime is None or regime == 'fsi':
		# 	print('starting FSI computation')
		# 	server = StructureServer('localhost', port, solver, pOut)  # , pRef, rhoRef, lRef, saveRate, animation)
		# elif regime == 'ale':
		# 	print('starting ALE computation')
		# 	server = StructureServerALE('localhost', port, solver)
		# elif regime == 'stationary':
		# 	print('starting test computation')
		# 	server = StructureServerTest('localhost', port, solver)
		# else:
		# 	raise ValueError('parameter regime has a wrong value')

		server.run()
	finally:
		os.chdir(currentPath)


if __name__ == '__main__':
	run()
