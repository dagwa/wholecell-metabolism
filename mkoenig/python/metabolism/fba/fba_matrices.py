"""
-----------------------------------------------------------------------------
Creating Matrices from the original matlab dump.
These matrices are representing the actual FBA problem exactly. They are 
calculated once at Karr model initialization and than used for the 
FBA calculation with dynamical updated fluxBounds.

@author: Matthias Koenig
@date: 2015-03-13
-----------------------------------------------------------------------------

   Knowledge Base
   ===============
   M. genitalium metabolism was reconstructed from a variety of sources,
   including flux-balance analysis metabolic models of other bacterial species
   and the reaction kinetics database SABIO-RK, and was organized into 504
   reactions in the knowledge base. These reactions are loaded into this process
   by the initializeConstants method.

      Object       No.
      ==========   ===
      substrates   585
      enzymes      104
      reactions    504
        chemical   ?
        transport  ?

   Representation
   ===============
   The properties substrates and enzymes represent the counts of metabolites
   and metabolic enzymes.

   fbaReactionStoichiometryMatrix represents the stoichiometry and compartments
   of metabolites and biomass in each of the 641 chemical/transport reactions,
   exchange pseudoreactions, and biomass production pseudoreaction.

   fbaReactionCatalysisMatrix represents the enzyme which catalyzes each
   reaction. 
   fbaEnzymeBounds represents the foward and backward kcat of the
   catalyzing enzyme of each reaction. 
   fbaReactionBounds represents the maximal import and export rates of each
   metabolite. 
   fbaObjective indicates which reaction represents the biomass production pseudoreaction. 
   fbaRightHandSide is a vector of zeros representing the change 
   in concentration over time of each metabolite and biomass. 
   metabolismProduction is redundant with the
   biomass production reaction in fbaReactionStoichiometryMatrix.
   metabolismProduction is calculated by summing the metabolic demands of all
   the other processes over the entire cell cycle. The table below lists the
   units of several properties of this process.

      Property                       Units
      ===========================    ==============================
      fbaEnzymeBounds                molecules/enzyme/s
      fbaReactionBounds              molecules/(gram dry biomass)/s
      substrates                     molecules
      enzymes                        molecules
      stepSizeSec                    s
      lowerBounds                    reactions/s
      upperBounds                    reactions/s
      growth                         cell/s
      biomassComposition             molecules/cell
      metabolismProduction           molecules/cell
      chamberVolume                  L
      setValues                      molecules/chamber
      growthAssociatedMaintanence
"""

from pandas import DataFrame
from pandas import Series

##############################################################################
# constant definitions in Karr model
##############################################################################
stepSizeSec = 1  # defined in Process
realmax = 1e6  # maximal flux bound
compartmentIndexs_cytosol = 1  # defined in Metabolism
compartmentIndexs_extracellular = 2  # defined in Metabolism
compartmentIndexs_membrane = 3  # defined in Metabolism

# Load state data from matlab dump (here also the constant matrices
# and indices are defined)
import scipy.io
import numpy as np
import os
from metabolism_settings import DATA_DIR, RESULTS_DIR
matrix_dir = os.path.join(RESULTS_DIR, 'fba_matrices')
matlab_data = os.path.join(DATA_DIR, 'matlab_dumps/Process_Metabolism.mat')

state = scipy.io.loadmat(matlab_data)
for key, value in sorted(state.iteritems()):
    if isinstance(value, np.ndarray):
        print key, value.shape
    else:
        print key

###############################################
# FULL SET
###############################################
# 645 Reactions
# 585 Substrates
# 104 Enzymes
# 3 Compartments
# reactionStoichiometryMatrix (585, 645, 3)
###############################################
# Original identifier of the matrices used in matlab.
# These matrices have to be broken down to 2D FBA matrix
rids = state["reactionWholeCellModelIDs"]
sids = state["substrateWholeCellModelIDs"]
eids = state["enzymeWholeCellModelIDs"]


def clean_whole_cell_Ids(ids, col_id):
    df = DataFrame(ids, columns=[col_id])
    for k in df.index:
        df.ix[k] = df.ix[k][0]
    df = df.set_index(df[col_id])
    return df

r_df = clean_whole_cell_Ids(rids, 'rid')
s_df = clean_whole_cell_Ids(sids, 'sid')
e_df = clean_whole_cell_Ids(eids, 'eid')

###############################################
# FBA SET
###############################################
# The full set of reactions and compounds is reduced to the information 
# needed for the FBA calculation. 
# Main challenge is to match the identifiers from the original larger set
# to the subset used in the FBA.
#
# 504 Reactions
# 376 Substrates
# 104 Enzymes


##########################################################################
# <Enzymes>
##########################################################################
# Find recursively the gene associations for given proteins.
# Every ProteinMonomer has an associated gene. Every ProteinComplex
# consists of ProteinMonomers and/or ProteinComplexes. 
from public.models import Protein, ProteinMonomer, ProteinComplex

def get_genes_for_protein(p):
    ''' Recursive finding of all genes for ProteinComplexes '''
    genes = set()
    if p.model_type == 'ProteinMonomer':
        pm = ProteinMonomer.objects.get(wid=p.wid)
        genes.add(pm.gene.wid)
        # print pm, pm.model_type, '=>', pm.gene.wid
        
    elif p.model_type == 'ProteinComplex':
        pc = ProteinComplex.objects.get(wid=p.wid)
        # find all the proteinmonomers involved
        for participant in pc.biosynthesis.all():
            if participant.molecule_id == pc.id:
                # do not further follow the protein we are already looking at
                continue
            else:
                protein = Protein.objects.get(id=participant.molecule_id)
                genes = genes.union(get_genes_for_protein(protein))
    return genes
        
# find genes
all_genes = []
for eid in e_df.eid:
    p = Protein.objects.get(wid=eid)
    genes = get_genes_for_protein(p)     
    print p, '|', p.model_type, '|', p.name
    print '=>', ','.join(genes)
    all_genes.append(','.join(genes)) # string format for DataFrame

e_df['genes'] = Series(all_genes, index=e_df.index)
e_df

##########################################################################
# <Reactions>
##########################################################################
# The reaction state vector consists of multiple parts, summing up to total
# 504 reaction elements:
# 336 Conversion
#   +  124 ExternalExchange
#   +   42 InternalExchange
#   +    1 biomassProduction
#   +    1 biomassExchange
#   =  504
# Which elements in the reaction vector are what type is defined in the fbaReactionIndexs matrices
# TODO: the naming of the external and internal exchange reactions should be adapted
# to the ones used in the provided Karr SBML model.
fbaReactionIndexs_metabolicConversion = state['fbaReactionIndexs_metabolicConversion']  # [336x1]
fbaReactionIndexs_metaboliteExternalExchange = state['fbaReactionIndexs_metaboliteExternalExchange']  # [124x1]
fbaReactionIndexs_metaboliteInternalExchange = state['fbaReactionIndexs_metaboliteInternalExchange']  # [42x1]
fbaReactionIndexs_biomassProduction = state['fbaReactionIndexs_biomassProduction']  # [1x1] (index 503)
fbaReactionIndexs_biomassExchange = state['fbaReactionIndexs_biomassExchange']  # [1x1] (index 504)

# The mapping of the 336 reaction ids from the large reaction vector to the fba subset
reactionIndexs_fba = state["reactionIndexs_fba"]  # [336,1]

# Reaction identifier
r_fba_df = DataFrame(range(0, 504), columns=['rid'])
# metabolic conversion
for k, item in enumerate(reactionIndexs_fba):
    index = item[0] - 1
    r_fba_df.loc[k] = r_df['rid'][index]

# metabolite external Exchange
for k, index in enumerate(fbaReactionIndexs_metaboliteExternalExchange):
    index -= 1
    r_fba_df.loc[index] = 'metaboliteExternalExchange_ex_{}'.format(k)

# metabolite internal Exchange
for k, index in enumerate(fbaReactionIndexs_metaboliteInternalExchange):
    index -= 1
    r_fba_df.loc[index] = 'metaboliteInternalExchange_ix_{}'.format(k)

# biomass production
index = fbaReactionIndexs_biomassProduction[0] - 1
r_fba_df.loc[index] = 'biomassProduction'.format(k)

# biomass consumption
index = fbaReactionIndexs_biomassExchange[0] - 1
r_fba_df.loc[index] = 'biomassConsumption'.format(k)

# set index
r_fba_df = r_fba_df.set_index(r_fba_df['rid'])

# Add reaction bound and objective information to the reaction vector
# maximal import and export rates (upper, lower)
fbaReactionBounds = state['fbaReactionBounds']  # [504x2]
r_fba_df['lb_fbaReactionBounds'] = fbaReactionBounds[:, 0]
r_fba_df['ub_fbaReactionBounds'] = fbaReactionBounds[:, 1]

# enzyme kinetics kcat (upper, lower)
fbaEnzymeBounds = state['fbaEnzymeBounds']  # [504x2]
np.isnan(fbaEnzymeBounds).any()
fbaEnzymeBounds[np.isnan(fbaEnzymeBounds[:, 0]), 0] = -np.inf
fbaEnzymeBounds[np.isnan(fbaEnzymeBounds[:, 1]), 1] = np.inf
np.isnan(fbaEnzymeBounds).any()

r_fba_df['lb_fbaEnzymeBounds'] = fbaReactionBounds[:, 0]
r_fba_df['ub_fbaEnzymeBounds'] = fbaReactionBounds[:, 1]

# fbaObjective indicates which reaction represents the biomass production pseudo-reaction.
fbaObjective = state['fbaObjective']  # [504x1]
r_fba_df['fbaObjective'] = fbaObjective

# save the reaction information
r_fba_df.to_csv(os.path.join(matrix_dir, 'r_fba.csv'), sep="\t", index=False)

##########################################################################
# <Substrates>
##########################################################################
# Calculation of the substrate identifier vector. Similar strategy than
# in the reaction id vector. Necessary to get the actual identifiers of the 
# FBA subset.
#
#   368 fba substrate
# +   7 internal exchange
# +   1 biomass
# sum = 376
fbaSubstrateIndexs_substrates = state['fbaSubstrateIndexs_substrates']  # [368x1]
fbaSubstrateIndexs_metaboliteInternalExchangeConstraints = state[
    'fbaSubstrateIndexs_metaboliteInternalExchangeConstraints']          # [7x1]
fbaSubstrateIndexs_biomass = state['fbaSubstrateIndexs_biomass']         # [1x1]
substrateIndexs_fba = state["substrateIndexs_fba"]                       # [376x1]

# substrate identifier
s_fba_df = DataFrame(range(0, 376), columns=['sid'])

# substrates
N_substrates = 585
index_offset = 8
for k, item in enumerate(substrateIndexs_fba):
    index = item[0] - 1
    compartment = 'c'
    if index >= N_substrates:
        compartment = 'e'
    if index >= 2*N_substrates:
        compartment = 'm'
    # get the columns via modulo
    index %= N_substrates

    # add the correct compartment    
    if compartment != 'c':
        sid = '{}__{}'.format(s_df['sid'][index], compartment)
    else:
        sid = '{}'.format(s_df['sid'][index])

    if compartment != 'c':
        # 8 elements (internal exchange and biomass) are out of order
        # so offset is nedessary
        k += index_offset
    s_fba_df.loc[k] = sid

# internal exchange
for k, index in enumerate(fbaSubstrateIndexs_metaboliteInternalExchangeConstraints):
    index -= 1
    s_fba_df.loc[index] = 'metaboliteInternalExchangeConstraints_ix_{}'.format(k)

# biomass
for k, index in enumerate(fbaSubstrateIndexs_biomass):
    index -= 1
    s_fba_df.loc[index] = 'biomass'.format(k)

# fbaRightHandSide is a vector of zeros representing the change 
# in concentration over time of each metabolite and biomass. 
fbaRightHandSide = state['fbaRightHandSide']  # [376x1]
s_fba_df['fbaRightHandSide'] = fbaRightHandSide

# set index and save
s_fba_df = s_fba_df.set_index(s_fba_df['sid'])
s_fba_df.to_csv(os.path.join(matrix_dir, 's_fba.csv'), sep="\t", index=False)

##########################################################################
# Reaction Stoichiometry Matrxix [376x504]
##########################################################################
# fbaReactionStoichiometryMatrix represents the stoichiometry and compartments
#   of metabolites and biomass in each of the 641 chemical/transport reactions,
#   exchange pseudoreactions, and biomass production pseudoreaction.
fbaReactionStoichiometryMatrix = state['fbaReactionStoichiometryMatrix']  # [376x504]

mat_stoichiometry = DataFrame(fbaReactionStoichiometryMatrix, columns=r_fba_df.index)
mat_stoichiometry = mat_stoichiometry.set_index(s_fba_df.index)
mat_stoichiometry.to_csv(os.path.join(matrix_dir, 'fbaReactionStoichiometryMatrix.csv'), sep="\t", index=False)

##########################################################################
# fbaReactionCatalysisMatrix  [504x104]
##########################################################################
# fbaReactionCatalysisMatrix represents the enzyme which catalyzes each
# reaction. 
fbaReactionCatalysisMatrix = state['fbaReactionCatalysisMatrix']  # [504x104]

mat_catalysis = DataFrame(fbaReactionCatalysisMatrix, columns=e_df.index)
mat_catalysis = mat_catalysis.set_index(r_fba_df.index)
mat_catalysis.to_csv(os.path.join(matrix_dir, 'fbaReactionCatalysisMatrix.csv'), sep="\t", index=False)
