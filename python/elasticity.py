import os
import sys
import subprocess

structureSolverPath = 'D:\Soubory\Ostatni\AlesTurek\elasticity'
sys.path.insert(1, structureSolverPath)

from dynamics import Dynamics
import statics
from parameters import Parameters, ParameterException, Model


paramFileName = 'structuralParameters.txt'

simulationDirPath = os.path.join(structureSolverPath, "simulations")

__geo__ = ''
__sim__ = 'default'


def loadArgs():
	"""Load geometry and simulation name from 'args.txt' file.
	loadArgs() returns two strings with the geometry and simulation name."""

	# file = open(os.path.join(structureSolverPath, 'args.txt'), 'r')
	# line = file.read()
	# line = line.strip()
	# strArray = line.split(" ")
	# return strArray[0], strArray[1]

	return __geo__, __sim__


def getPath(*argv):
	"""Return mesh, simulation and output paths.
	[meshPath, simulPath, outputPath] = getPath returns path of the mesh,
	simulation and outpud folders for the current simulation.
	[meshPath, simulPath, outputPath] = getPath(geometry, simulation) returns
	path of the mesh, simulation nad outpud folders for the simulation
	specified by the name of the geometry and simulation."""

	if len(argv) == 0:
		geometry, simulation = loadArgs()

	if len(argv) == 1:
		geometry, simulation = loadArgs()
		simulation = "default"

	meshPath = os.path.join(structureSolverPath, 'simulations', geometry, 'mesh')
	simulPath = os.path.join(structureSolverPath, 'simulations', geometry, simulation)
	outputPath = os.path.join(simulPath, 'output')
	animPath = os.path.join(simulPath, 'animation')

	if not os.path.isdir(meshPath):
		print("Geometry '" + geometry + "' does not exist.")
		return

	if not os.path.isdir(simulPath):
		print("Geometry '%s' does not contain simulation '%s'" % (geometry, simulation))
		return

	# if not os.path.isdir(outputPath):
	# 	outputPath = ""

	return meshPath, simulPath, outputPath, animPath


def getParam():
	geoPath, simPath, outPath = getPath()
	paramFilePath = os.path.join(simPath, paramFileName)

	with open(paramFilePath) as paramFile:
		lines = paramFile.read().splitlines()

	paramDct = {}
	for line in lines:
		line = line.split('%')[0]
		line = line.split('#')[0]
		tokens = line.split('=')
		if len(tokens) == 2:
			command = tokens[0].strip()
			value = tokens[1].strip()
			paramDct[command] = value

	return paramDct


def args(geometry=None, simulation='default'):
	"""Change or print geometry and simulation name.

	args(geometry, simulation='default') changes current geometry and simulation.
	It does this by saving the geometry and simulation names into 'args.txt'.
	When FlowPro is executed, it will run this simulation.

	args(geometry) does the above and sets simulation name as 'default'.

	args() prints current geometry and simulation name."""

	if geometry is None:
		geometry, simulation = loadArgs()
		print("Geometry name:", geometry)
		print("Simulation name:", simulation)
		return

	geometry = os.path.relpath(geometry)

	geomPath = os.path.join(structureSolverPath, 'simulations', geometry)
	if not os.path.isdir(geomPath):
		print("Geometry '" + geometry + "' does not exist.")
		return

	simPath = os.path.join(geomPath, simulation)
	if not os.path.isdir(simPath):
		print("Geometry '%s' does not contain simulation '%s'" % (geometry, simulation))
		print("Simulation will be created after running command 'param'")

	# file = open(os.path.join(structureSolverPath, 'args.txt'), 'w')
	# file.write(geometry + ' ' + simulation)
	# file.close()

	global __geo__, __sim__
	__geo__ = geometry
	__sim__ = simulation


def listSims():
	"""Print the list of geometries and simulations."""

	__listSubdir('')


def __listSubdir(geoName):
	geoPath = os.path.join(simulationDirPath, geoName)
	subdirs = [f for f in os.listdir(geoPath) if os.path.isdir(os.path.join(geoPath, f))]
	if 'mesh' in subdirs:
		subdirs.remove('mesh')
		print(geoName + ': ' + ', '.join(subdirs))
	else:
		for subdir in subdirs:
			__listSubdir(os.path.join(geoName, subdir))


def run():
	"""Run Structure Solver."""
	geoPath, simPath, outPath, animation = getPath()
	parameters = Parameters(simPath)

	if parameters.model == Model.statics:
		statics.solve(geoPath, simPath, outPath)
	elif parameters.model == Model.dynamics:
		solver = Dynamics(geoPath, simPath, outPath, animation)
		solver.solve()


if __name__ == '__main__':
	args('turekhron/2', 'csm3')
	dct = getParam()
	print(dct)


def info():
	"""List all available methods along with the first line from their docs."""
	import inspect
	import elasticity as ela

	print(inspect.getdoc(ela))
	print()
	print("Available methods:")
	print()

	publicMethods = [s for s in dir(ela) if callable(getattr(ela, s))]
	for method in publicMethods:
		doc = inspect.getdoc(getattr(ela, method))
		if doc is None:
			doc = ''
		firstLine = doc.split('\n')[0]
		print("%-12s  %s" % (method, firstLine))


def param():
	"""Open text file with parameters for the problem."""

	geometry, simulation = loadArgs()
	geomPath = os.path.join(structureSolverPath, 'simulations', geometry)

	if not os.path.isdir(geomPath):
		print("geometry does not exist")
		return

	paramPath = os.path.join(geomPath, simulation, paramFileName)

	if not os.path.isfile(paramPath):
		answer = input("Simulation " + simulation + " does not exist, do you want to creat it [Y/n]? ")

		ans = answer.lower()
		if ans == "y" or ans == "yes" or ans == "":
			subprocess.call('java -jar FlowProManager.jar createparamfile', shell=True, cwd=structureSolverPath)
		else:
			return

	os.system("subl " + paramPath)
