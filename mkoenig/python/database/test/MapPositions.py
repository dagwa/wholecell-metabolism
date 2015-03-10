'''
Get all the position information for the metabolic map from the 
MapPositions in the Knowledge Base.

@author: Matthias Koenig
@date: 2015-03-09
'''

from pandas import DataFrame
import pandas as pd

from public.models import Metabolite, Reaction
from public.models import MetaboliteMapCoordinate, ReactionMapCoordinate

#------------------------
# Metabolite positions
#------------------------
m_df = DataFrame(columns=('id', 'model_type', 'wid', 'name', 'compartment', 'x', 'y'))
m_cor = MetaboliteMapCoordinate.objects.all()
for k, p in enumerate(m_cor):
    # get the associated metabolite (should be exactly 1) 
    m = p.metabolites.all()[0]
    # add metabolite and position
    m_df.loc[k] = (m.id, m.model_type, m.wid, m.name,
                   p.compartment, p.x, p.y)
m_df = m_df.set_index(m_df.id)
# set the data types
m_df[['id', 'x', 'y']] = m_df[['id', 'x', 'y']].astype(int)

m_df.head(10)


#------------------------
# Reaction positions
#------------------------
r_df = DataFrame(columns=('id', 'model_type', 'wid', 'name', 
                          'path', 'value_x', 'value_y', 'label_x', 'label_y'))
r_cor = ReactionMapCoordinate.objects.all()
for k, p in enumerate(r_cor):
    # get reaction
    r = p.reactions.all()[0]
    r_df.loc[k] = (r.id, r.model_type, r.wid, r.name, 
                   p.path, p.value_x, p.value_y, p.label_x, p.label_y)
r_df = r_df.set_index(r_df.id)
# set the data types
r_df[['id', 'value_x', 'value_y', 'label_x', 'label_y']] = r_df[['id', 'value_x', 'value_y', 'label_x', 'label_y']].astype(int)

r_df.head(10)
