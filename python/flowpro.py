"""Python interface for FlowPro."""

import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt


flowProPath = 'E:\\workspace\\FlowProPackage\\FlowPro'
paramFileName = 'parameters.txt'

simulationDirPath = os.path.join(flowProPath, "simulations")


def loadArgs():
	"""Load geometry and simulation name from 'args.txt' file.
	loadArgs() returns two strings with the geometry and simulation name."""

	file = open(os.path.join(flowProPath, 'args.txt'), 'r')
	line = file.read()
	line = line.strip()
	strArray = line.split(" ")
	return strArray[0], strArray[1]


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

	meshPath = os.path.join(flowProPath, 'simulations', geometry, 'mesh')
	simulPath = os.path.join(flowProPath, 'simulations', geometry, simulation)
	outputPath = os.path.join(simulPath, 'output')

	if not os.path.isdir(meshPath):
		print("Geometry " + geometry + " does not exist.")
		return

	if not os.path.isdir(simulPath):
		print("Simulation " + geometry + " does not exist.")
		return

	if not os.path.isdir(outputPath):
		outputPath = ""

	return meshPath, simulPath, outputPath


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

	geomPath = os.path.join(flowProPath, 'simulations', geometry)
	if not os.path.isdir(geomPath):
		print("Geometry " + geometry + " does not exist.")
		return

	simPath = os.path.join(geomPath, simulation)
	if not os.path.isdir(simPath):
		print("Simulation " + geometry + " " + simulation + " does not exist")
		print("Simulation will be created after running command 'param'")

	file = open(os.path.join(flowProPath, 'args.txt'), 'w')
	file.write(geometry + ' ' + simulation)
	file.close()


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


def run(*argv):
	"""Run FlowPro."""

	if len(argv) == 0:
		subprocess.call(['start', 'run.bat'], shell=True, cwd=flowProPath)
	else:
		subprocess.call('java -d64 -Xmx8g -jar FlowPro.jar master' + str(argv[0]), shell=True, cwd=flowProPath)


def show(*argv):
	"""Export data for post-processing.
	show('<variableName1, variableName2, ...> -i<iteration> -f<format>')

	formats: vtk (Paraview), txt (plain text)

	variable names: depends on the model used

	example:

	  show('mach pressure -i100 -ftxt')"""

	command = "java -d64 -Xmx8g -jar FlowPro.jar postprocessing "
	for str in argv:
		command = command + " " + str

	subprocess.call(command, shell=True, cwd=flowProPath)


	try:
		paths = getPath()
		vertices = np.loadtxt(paths[2] + "vertices.txt")
		#if(size(vertices,2) > 2)
		#    return
		#end
		elements = np.loadtxt(paths[0] + "elements.txt")
		elementType = np.loadtxt(paths[0] + "elementType.txt")

		for arg in argv:
			if not arg[0] == " ":
				quantity = np.loadtxt(paths[2] + arg + ".txt")
				myContour(elements, vertices, quantity, arg)
		plt.show()

	except:
		print("Unable to show results in python!")


def showResiduum():
	paths = getPath()
	residuum = np.loadtxt(paths[1] + "residuum.txt")
	plt.figure()
	plt.semilogy(residuum[:,1],linewidth=2.0,color='k')
	plt.xlabel('resid')
	plt.ylabel('iteration')
	plt.grid()
	plt.show()


def fetcher(*argv):
	command = "java -jar Fetcher.jar"
	for str in argv:
		command = command + " " + str

	subprocess.call(command, shell=True, cwd=flowProPath)


def info():
	"""List all available methods along with the first line from their docs."""
	import inspect
	import flowpro as fp

	print(inspect.getdoc(fp))
	print()
	print("Available methods:")
	print()

	publicMethods = [s for s in dir(fp) if callable(getattr(fp, s))]
	for method in publicMethods:
		doc = inspect.getdoc(getattr(fp, method))
		if doc is None:
			doc = ''
		firstLine = doc.split('\n')[0]
		print("%-12s  %s" % (method, firstLine))


def myContour(triangles, vertices, quantity, name):
	#plt.triplot(vertices[:,0], vertices[:,1], triangles)
	plt.figure()
	plt.tricontourf(vertices[:,0], vertices[:,1], triangles, quantity, 20)
	plt.gca().set_aspect('equal')
	plt.colorbar()
	plt.title(name)
	plt.xlabel('x')
	plt.ylabel('y')


def param():
	"""Open text file with parameters for the problem."""

	geometry, simulation = loadArgs()
	geomPath = os.path.join(flowProPath, 'simulations', geometry)

	if not os.path.isdir(geomPath):
		print("geometry does not exist")
		return
		# raise Exception("geometry does not exist")

	paramPath = os.path.join(geomPath, simulation, "parameters.txt")

	if not os.path.isfile(paramPath):
		answer = input("Simulation " + simulation + " does not exist, do you want to creat it [Y/n]? ")

		ans = answer.lower()
		if ans == "y" or ans == "yes" or ans == "":
			subprocess.call('java -jar FlowProManager.jar createparamfile', shell=True, cwd=flowProPath)
		else:
			return

	os.system("subl " + paramPath)
