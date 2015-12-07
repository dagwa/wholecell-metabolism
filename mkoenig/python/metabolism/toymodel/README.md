# Toy model
Simple toy model to demonstrate the coupling of SBML models with different modelling frameworks.

## Modelling framework
For the simulator to understand with which modelling framework submodels/models should be simulated it is
necessary to annotate the `<model>` element with the modelling framework. For the annotation the SBO terms below the
`SBO:0000004 - modelling framework` are used.

```
SBO:0000004 - modelling framework
|    SBO:0000062 continuous framework
|    SBO:0000063 discrete framework
|    SBO:0000624 flux balance framework
|    SBO:0000234 logical framework
```
In a first version the SBOTerm of the model is set. This is currently not legal, but the simplest approach.

## Submodels
![submodel overview](docs/toymodel_overview.png)



## Create models and simulate
The following content is available for the toy model
```
toymodel
|   settings.py: file names and destinations
|   model_factory.py: creates the individual submodels
|   comp_factory.py: creates the combined comp models
|   simulator.py: simulator for the comp model
|   run_all.py: script
|   models/: created models
|   docs/: documentation of models (reports)
|   results/: results of simulations
```

### Requirements
Build the following libraries with their python bindings from the latest source
* `libsbml`
* `roadrunner`
* `cobrapy`
* 'sbmlutils' - python sbml utils (package missing)

### Create models & simulate
```
python run_all.py
```



