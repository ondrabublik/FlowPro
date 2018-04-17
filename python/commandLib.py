import sys, os
import numpy as np
import matplotlib.pyplot as plt

def loadArgs():
    file = open("../args.txt", "r") 
    line = file.read()
    line = line.strip()
    strArray = line.split(" ")
    return strArray;

def getPath(*argv):
    if len(argv) == 0:
        strArray = loadArgs()
        geometry = strArray[0]
        simulation = strArray[1]
    
    if len(argv) == 1:
        strArray = loadArgs()
        geometry = strArray[0]
        simulation = "default"

    meshPath = "../simulations/" + geometry + "/mesh/"
    simulPath = "../simulations/" + geometry + "/" + simulation + "/"
    outputPath = "../simulations/" + geometry + "/" + simulation + "/output/"

    if not os.path.isdir(meshPath):
        print "Geometry " + geometry + " does not exist."
        return

    if not os.path.isdir(simulPath):
        print "Simulation " + geometry + " does not exist."
        return

    if not os.path.isdir(outputPath):
        outputPath = "";
    
    return [meshPath,simulPath,outputPath]


def args(*argv):
    if len(argv) == 0:
        strArray = loadArgs()
        geometry = strArray[0]
        simulation = strArray[1]
        print "Geometry name:", geometry
        print "Simulation name:", simulation
        return
    elif len(argv) == 1:
        print "Setting simulation to default."
        geometry = argv[0]
        simulation = "default"
    else:
        geometry = argv[0]
        simulation = argv[1]

    path = "../simulations/" + geometry + "/"
    if not os.path.isdir(path):
        print "Geometry " + geometry + " does not exist."
        return

    path = path + simulation
    if not os.path.isdir(path):
        print "Simulation " + geometry + " " + simulation + " does not exist";
        print "the simulation will be created after running the command 'param'"

    file = open("../args.txt","w") 
    file.write(geometry + " " + simulation)
    file.close()

def list():
    dirName = "../simulations";
    dir = os.listdir(dirName)
    for subdir in dir:
        subdirs = os.listdir(dirName + "/" + subdir)
        text = subdir + ":"
        for str in subdirs:
            text = text + " " + str
        print text

def run(*argv):
    if len(argv) == 0:
        os.chdir("..")
        os.system("start java -d64 -Xmx8g -jar FlowPro.jar local")
        os.chdir("python")
    else:
        os.chdir("..")
        os.system("java -d64 -Xmx8g -jar FlowPro.jar master " + str(argv[0]))
        os.chdir("python")

def show(*argv):
    command = "start java -d64 -Xmx8g -jar FlowPro.jar postprocessing "
    for str in argv:
        command = command + " " + str
    os.chdir("..")
    os.system(command)
    os.chdir("python")

    try:
        paths = getPath();
        vertices = np.loadtxt(paths[2] + "vertices.txt")
        #if(size(vertices,2) > 2)
        #    return;
        #end
        elements = np.loadtxt(paths[0] + "elements.txt")
        elementType = np.loadtxt(paths[0] + "elementType.txt")

        for arg in argv:
            if not arg[0] == " ":
                quantity = np.loadtxt(paths[2] + arg + ".txt")
                myContour(elements, vertices, quantity, arg);
        plt.show()

    except:
        print "Unable to show results in python!"
    
def showResiduum():
    paths = getPath();
    residuum = np.loadtxt(paths[1] + "residuum.txt")
    plt.figure()
    plt.semilogy(residuum[:,1],linewidth=2.0,color='k');
    plt.xlabel('resid')
    plt.ylabel('iteration')
    plt.grid()
    plt.show()
    

def fetcher(*argv):
    command = "java -jar Fetcher.jar"
    for str in argv:
        command = command + " " + str
    os.chdir("..")
    os.system(command)
    os.chdir("python")

def info():
    print "**********************************"
    print "*                                *"
    print "*    FlowPro python interface    *"
    print "*                                *"
    print "**********************************"
    print " "
    print "commands available: "
    print "args() - for switching actual simulation"
    print "list() - list of simulations"
    print "show() - list of simulations"
    print "fetcher() - "

def myContour(triangles, vertices, quantity, name):
    #plt.triplot(vertices[:,0], vertices[:,1], triangles)
    plt.figure()
    plt.tricontourf(vertices[:,0], vertices[:,1], triangles, quantity, 20)
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.title(name)
    plt.xlabel('x')
    plt.ylabel('y')












