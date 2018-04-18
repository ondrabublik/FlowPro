# FlowPro

FlowPro is a numerical tool for solving a general systems of hyperbolic partial differential equations
 based on discontinuous Galerkin finite element method (DGFEM). 
The main feature of the software is that you can easily define your own mathematical model.

## Getting Started


### Prerequisites
-java 8 is installed

## Running the tests
At first, you must set a path of the simulation in file args.txt Open file args.txt, write "examples/NACA default" and save it. 
Then run the comand at commandline

java -jar FlowPro.jar local

To show results run at commandline

java -jar FlowPro.jar postprocessing mach pressure -fvtk

## Built With

* [Ant] - build.xml

## Contributing

## Authors
Ondrej Bublik [obublik@kme.zcu.cz]
Ales Pecka [pecka@kme.zcu.cz]
Jan Vimmr [jvimmr@kme.zcu.cz]

## License

This project is licensed under the BSD License - see the [license.txt](license.txt) file for details

## Acknowledgments

