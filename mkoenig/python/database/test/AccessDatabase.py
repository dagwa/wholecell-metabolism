'''
Test the object mapping.
@author: mkoenig
@date: 2015-03-09
'''

from public.models import Metabolite, Reaction, Gene

# print some metabolites
metabolites = Metabolite.objects.all()[:10]    
for m in metabolites: 
    print m.pk, m.name, m.charge
# find some genes
Gene.objects.filter(name__contains="2")

# find some reactions
rs = Reaction.objects.filter(name__startswith='A')
print rs
# print the stoichiometry of first reaction
print rs[0].stoichiometry
print rs[0].stoichiometry.all()

from public.models import Entry, SpeciesComponent
from pandas import DataFrame
import pandas as pd

# write pandas DataFrame


entries = Entry.objects.all()[:10]
entries_df = DataFrame(columns=('id', 'model_type', 'wid', 'name'))
for k, e in enumerate(entries):
    entries_df.loc[k] = (e.id, e.model_type, e.wid, e.name)
entries_df = entries_df.set_index(entries_df.id)
print entries_df

scs = SpeciesComponent.objects.all()[:10]
print scs
print [s.parent_ptr_entry for s in scs]

# Create the positions for the layouts
from public.models import Metabolite, Reaction
from public.models import MetaboliteMapCoordinate, ReactionMapCoordinate

xys = MetaboliteMapCoordinate.objects.all()[:10]
for xy in xys:
    print xy.id, xy.compartment, xy.x, xy.y
    m = xy.metabolites.all()[0]
    print m
    print m.id, m.model_type, m.wid, m.name
    
xys = ReactionMapCoordinate.objects.all()[:10]
for xy in xys:
    print xy.id, xy.path, xy.value_x, xy.value_y, xy.label_x, xy.label_y
    r = xy.reactions.all()[0]
    print r
    print r.id, r.model_type, r.wid, r.name


