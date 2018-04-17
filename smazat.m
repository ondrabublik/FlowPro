function smazat
copyfile('../FlowProModules/Aerodynamics/dist/Aerodynamics.jar','modules/equations/Aerodynamics.jar')
copyfile('../FlowProModules/IdealMHD/dist/IdealMHD.jar','modules/equations/IdealMHD.jar')
copyfile('../FlowProModules/LabyrinthSeal3D/dist/LabyrinthSeal3D.jar','modules/equations/LabyrinthSeal3D.jar')
copyfile('../FlowProModules/Multiphase/dist/Multiphase.jar','modules/equations/Multiphase.jar')
copyfile('../FlowProModules/ScalarEquations/dist/ScalarEquations.jar','modules/equations/ScalarEquations.jar')
copyfile('../FlowProModules/ShallowWater/dist/ShallowWater.jar','modules/equations/ShallowWater.jar')
copyfile('../FlowProModules/RigidBody/dist/RigidBody.jar','modules/dynamics/RigidBody.jar')
copyfile('../FlowProManager/dist/FlowProManager.jar','FlowProManager.jar')
copyfile('../FlowProModules/SolutionMonitor/dist/SolutionMonitor.jar','modules/solutionmonitor/SolutionMonitor.jar')
% copyfile('../FlowProModules/ScalarEquations/dist/ScalarEquations.jar','modules/equations/ScalarEquations.jar')
% copyfile('../FlowProModules/Multiphase/dist/Multiphase.jar','modules/equations/Multiphase.jar')
!java -jar FlowProManager.jar update