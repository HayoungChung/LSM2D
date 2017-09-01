# LSM2D

## introduction
to get cythonized library of level-set toolbox, 

```bash
./lib/lsm2d_classwise
sh ./run_setup (case of bash shell)
cp ./build/lib*/lsm_classes.so ./../../LSTO2D
```

```iPython
run run_LSTO.run
```

## Error message:
``` Python
/home/hac210/Dropbox/packages/topOpt_MDO/LSTO2D/lsm2d_SLP_Group.pyc in setup(self)
     48         # constraint setup
     49         comp = ConstraintComp(lsm_solver = lsm_solver,
---> 50                              num_bpts = num_bpts)
     51         self.add_subsystem('constraint_comp', comp)
     52         comp.add_constraint('constraint')
"Key 'lsm_solver' cannot be set because it has not been declared."
```
>what puzzles me: these components share same commands, but only ConstraintComp() goes into trouble.

``` Python
    # distance computation
    comp = DistanceComp(lsm_solver = lsm_solver, 
                        num_dvs = num_dvs, 
                        num_bpts = num_bpts)
    self.add_subsystem('distance_comp',comp)
    self.connect('distance_comp.distances', 'constraint_comp.distances')
    self.connect('distance_comp.distances', 'objective_comp.distances')

    # objective setup
    comp = ObjectiveComp(lsm_solver = lsm_solver,
                            num_bpts = num_bpts)
    self.add_subsystem('objective_comp', comp)
    comp.add_objective('objective')

            
    # constraint setup
    comp = ConstraintComp(lsm_solver = lsm_solver,
                            num_bpts = num_bpts)
    self.add_subsystem('constraint_comp', comp)
    comp.add_constraint('constraint')
```

## MSG. to John:
this version of LSTO is the simplest version among I've tried: 
in this version, only 3 components are exposed that compute: 

    1. displacement
        - input: lambdas
        - outputs: displacement
    2. objective function
        - input: displacement
        - output: objective
    3.  constraint function
        - input: displacement
        - output: constraint

FYI, a legacy of different versions can be found in ./visualize folder :)

