'''
Test the object mapping.
@author: mkoenig
@date: 2015-03-09
'''

# import the objects you are interested in from models
# Detailed DB schema in schema & relationships in models.py
from public.models import Metabolite

# get first 10 metabolites
metabolites = Metabolite.objects.all()[:10]    
for m in metabolites: 
    # Metabolite inherits from Entry, so all fields are available
    print m.pk, m.id, m.model_type, m.wid, m.name, m.charge
    
# find all genes which have a 2 in it
from public.models import Gene
genes = Gene.objects.filter(name__contains="2")
print genes

# find all reactions starting with A
from public.models import Reaction
rs = Reaction.objects.filter(name__startswith='A')
print rs
# print the stoichiometry of first reaction
print rs[0].stoichiometry.all()

from public.models import Entry
from pandas import DataFrame
import pandas as pd

# create pandas DataFrame for all entries
entries = Entry.objects.all()
entries_df = DataFrame(columns=('id', 'model_type', 'wid', 'name'))
for k, e in enumerate(entries):
    # Bad in place extension of DataFrame (do not do in production code)
    entries_df.loc[k] = (e.id, e.model_type, e.wid, e.name)
entries_df = entries_df.set_index(entries_df.id)
# print the first 10
print entries_df.head(10)