# FlowPro

FlowPro is a numerical tool for solving a general systems of hyperbolic partial differential equations
 based on discontinuous Galerkin finite element method (DGFEM). 
The main feature of the software is that you can easily define your own mathematical model.

## Getting Started
1) download or clone FlowPro
2) compile FlowPro (run ant)
3) download or clone FlowProManager
4) compile FlowProManager (run ant)
5) copy FlowProManager.jar to same place as FlowPro.jar
6) run in terminal: java -jar FlowProManager.jar update (this update a manifest file in FlowPro.jar)

### Prerequisites
* java 8
* paraview (for visualisation) 

## Running the tests
At first, you must set a path of the simulation in file args.txt Open file args.txt, write "examples/NACA default" and save it. 
Then run the command at commandline

java -jar FlowPro.jar local

To show results run at commandline

java -jar FlowPro.jar postprocessing mach pressure -fvtk

## Built With

* [Ant] - build.xml

## Authors
Ondrej Bublik [obublik@kme.zcu.cz]
Ales Pecka [pecka@kme.zcu.cz]
Jan Vimmr [jvimmr@kme.zcu.cz]

## License

This project is licensed under the BSD License - see the [license.txt](license.txt) file for details

## Acknowledgments

