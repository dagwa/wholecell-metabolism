"""
Test interaction with gene associations.
"""

from __future__ import print_function, division
import libsbml


def read_gene_associations_from_fbc(sbml):
    """ Reading all the GeneAssociations from FBC. """
    doc = libsbml.readSBML(sbml)
    
    if doc.getPlugin("fbc") is not None:
        model_sbml = doc.getModel()
        mplugin = model_sbml.getPlugin("fbc");
        for ga in mplugin.getListOfGeneAssociations():            
             
            ass = ga.getAssociation()
            infix = ass.toInfix()  # get the rule string for cobray
            print(infix, ':', ass.toSBML())


##############################################################################
if __name__ == "__main__":
    import os
    sbml = os.path.join('/home/mkoenig/wholecell-metabolism/mkoenig/results', "Metabolism_annotated_4_L3V1.xml")
    read_gene_associations_from_fbc(sbml)
