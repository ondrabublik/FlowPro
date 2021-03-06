# FlowPro

FlowPro is a numerical tool for solving a general systems of hyperbolic partial differential equations
 based on discontinuous Galerkin finite element method (DGFEM). 
The main feature of the software is that you can easily define your own mathematical model.

## Getting Started
1) download or clone [FlowProPackage](https://github.com/ondrabublik/FlowProPackage)
2) compile FlowPro (gradle)

### Prerequisites
* java 8
* paraview (for visualisation) 

## Running first simulation
Open command line window and navigate to your FlowPro directory. FlowPro is executed with the command

java -jar FlowPro.jar master 0

When the simulation is done, you can generate results by running

java -jar FlowPro.jar postprocessing mach pressure -fvtk

In order to display the results, open file <FlowProDirectory>/simulations/examples/GAMM/default/output/results.vtk in paraview.

## Built With
* [Gradle](https://gradle.org/)

## Authors
Ondrej Bublik [obublik@kme.zcu.cz]  
Ales Pecka [pecka@kme.zcu.cz]  
Jan Vimmr [jvimmr@kme.zcu.cz]

## License

This project is licensed under the BSD License - see the [license.txt](license.txt) file for details

## Acknowledgments

