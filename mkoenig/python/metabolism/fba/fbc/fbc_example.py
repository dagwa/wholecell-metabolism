"""
Creates a simple fbc example network.
"""

from libsbml import *
from sbmlutils.validation import validate_sbml


def create_fbc_example():
    """ Create FBC example network.

    :return:
    :rtype:
    """
    sbmlns = SBMLNamespaces(3, 1, "fbc", 1)

    # create the document

    document = SBMLDocument(sbmlns)
    document.setPackageRequired("fbc", False)
    # create the Model
    model = document.createModel()

    # create the Compartment

    compartment = model.createCompartment()
    compartment.setId("compartment")
    compartment.setConstant(True)
    compartment.setSize(1)

    # create the Species

    species = model.createSpecies()
    species.setId("Node1")
    species.setCompartment("compartment")
    species.setBoundaryCondition(False)
    species.setConstant(False)
    species.setHasOnlySubstanceUnits(False)

    species = model.createSpecies()
    species.setId("Node2")
    species.setCompartment("compartment")
    species.setBoundaryCondition(False)
    species.setConstant(False)
    species.setHasOnlySubstanceUnits(False)

    species = model.createSpecies()
    species.setId("Node3")
    species.setCompartment("compartment")
    species.setBoundaryCondition(False)
    species.setConstant(False)
    species.setHasOnlySubstanceUnits(False)

    species = model.createSpecies()
    species.setId("Node4")
    species.setCompartment("compartment")
    species.setBoundaryCondition(False)
    species.setConstant(False)
    species.setHasOnlySubstanceUnits(False)

    species = model.createSpecies()
    species.setId("Node5")
    species.setCompartment("compartment")
    species.setBoundaryCondition(False)
    species.setConstant(False)
    species.setHasOnlySubstanceUnits(False)

    species = model.createSpecies()
    species.setId("Node6")
    species.setCompartment("compartment")
    species.setBoundaryCondition(False)
    species.setConstant(False)
    species.setHasOnlySubstanceUnits(False)

    species = model.createSpecies()
    species.setId("Node7")
    species.setCompartment("compartment")
    species.setBoundaryCondition(False)
    species.setConstant(False)
    species.setHasOnlySubstanceUnits(False)

    species = model.createSpecies()
    species.setId("Node8")
    species.setCompartment("compartment")
    species.setBoundaryCondition(False)
    species.setConstant(False)
    species.setHasOnlySubstanceUnits(False)

    species = model.createSpecies()
    species.setId("Node0")
    species.setCompartment("compartment")
    species.setBoundaryCondition(True)
    species.setConstant(False)
    species.setHasOnlySubstanceUnits(False)

    species = model.createSpecies()
    species.setId("Node9")
    species.setCompartment("compartment")
    species.setBoundaryCondition(True)
    species.setConstant(False)
    species.setHasOnlySubstanceUnits(False)

    # create reactions
    reaction = model.createReaction()
    reaction.setId("J0")
    reaction.setReversible(False)
    reaction.setFast(False)
    reactant = reaction.createReactant()
    reactant.setSpecies("Node0")
    reactant.setStoichiometry(1)
    reactant.setConstant(True)
    product = reaction.createProduct()
    product.setSpecies("Node1")
    product.setStoichiometry(1)
    product.setConstant(True)

    reaction = model.createReaction()
    reaction.setId("J1")
    reaction.setReversible(False)
    reaction.setFast(False)
    reactant = reaction.createReactant()
    reactant.setSpecies("Node1")
    reactant.setStoichiometry(1)
    reactant.setConstant(True)
    product = reaction.createProduct()
    product.setSpecies("Node2")
    product.setStoichiometry(1)
    product.setConstant(True)

    reaction = model.createReaction()
    reaction.setId("J2")
    reaction.setReversible(False)
    reaction.setFast(False)
    reactant = reaction.createReactant()
    reactant.setSpecies("Node2")
    reactant.setStoichiometry(1)
    reactant.setConstant(True)
    product = reaction.createProduct()
    product.setSpecies("Node3")
    product.setStoichiometry(1)
    product.setConstant(True)

    reaction = model.createReaction()
    reaction.setId("J3")
    reaction.setReversible(False)
    reaction.setFast(False)
    reactant = reaction.createReactant()
    reactant.setSpecies("Node1")
    reactant.setStoichiometry(1)
    reactant.setConstant(True)
    product = reaction.createProduct()
    product.setSpecies("Node4")
    product.setStoichiometry(1)
    product.setConstant(True)

    reaction = model.createReaction()
    reaction.setId("J4")
    reaction.setReversible(False)
    reaction.setFast(False)
    reactant = reaction.createReactant()
    reactant.setSpecies("Node4")
    reactant.setStoichiometry(1)
    reactant.setConstant(True)
    product = reaction.createProduct()
    product.setSpecies("Node3")
    product.setStoichiometry(1)
    product.setConstant(True)

    reaction = model.createReaction()
    reaction.setId("J5")
    reaction.setReversible(False)
    reaction.setFast(False)
    reactant = reaction.createReactant()
    reactant.setSpecies("Node3")
    reactant.setStoichiometry(1)
    reactant.setConstant(True)
    product = reaction.createProduct()
    product.setSpecies("Node5")
    product.setStoichiometry(1)
    product.setConstant(True)

    reaction = model.createReaction()
    reaction.setId("J6")
    reaction.setReversible(False)
    reaction.setFast(False)
    reactant = reaction.createReactant()
    reactant.setSpecies("Node5")
    reactant.setStoichiometry(1)
    reactant.setConstant(True)
    product = reaction.createProduct()
    product.setSpecies("Node6")
    product.setStoichiometry(1)
    product.setConstant(True)

    reaction = model.createReaction()
    reaction.setId("J7")
    reaction.setReversible(False)
    reaction.setFast(False)
    reactant = reaction.createReactant()
    reactant.setSpecies("Node6")
    reactant.setStoichiometry(1)
    reactant.setConstant(True)
    product = reaction.createProduct()
    product.setSpecies("Node7")
    product.setStoichiometry(1)
    product.setConstant(True)

    reaction = model.createReaction()
    reaction.setId("J8")
    reaction.setReversible(False)
    reaction.setFast(False)
    reactant = reaction.createReactant()
    reactant.setSpecies("Node5")
    reactant.setStoichiometry(1)
    reactant.setConstant(True)
    product = reaction.createProduct()
    product.setSpecies("Node8")
    product.setStoichiometry(1)
    product.setConstant(True)

    reaction = model.createReaction()
    reaction.setId("J9")
    reaction.setReversible(False)
    reaction.setFast(False)
    reactant = reaction.createReactant()
    reactant.setSpecies("Node8")
    reactant.setStoichiometry(1)
    reactant.setConstant(True)
    product = reaction.createProduct()
    product.setSpecies("Node7")
    product.setStoichiometry(1)
    product.setConstant(True)

    reaction = model.createReaction()
    reaction.setId("J10")
    reaction.setReversible(False)
    reaction.setFast(False)
    reactant = reaction.createReactant()
    reactant.setSpecies("Node7")
    reactant.setStoichiometry(1)
    reactant.setConstant(True)
    product = reaction.createProduct()
    product.setSpecies("Node9")
    product.setStoichiometry(1)
    product.setConstant(True)

    # add FBC constraints
    mplugin = model.getPlugin("fbc")

    bound = mplugin.createFluxBound()
    bound.setId("bound1")
    bound.setReaction("J0")
    bound.setOperation("equal")
    bound.setValue(10)

    objective = mplugin.createObjective()
    objective.setId("obj1")
    objective.setType("maximize")

    mplugin.setActiveObjectiveId("obj1")

    fluxObjective = objective.createFluxObjective()
    fluxObjective.setReaction("J8")
    fluxObjective.setCoefficient(1)

    sbml_out = "fbc_example.xml"
    writeSBMLToFile(document, sbml_out)

    validate_sbml(sbml_out, ucheck=False)

if __name__ == "__main__":
    create_fbc_example()
