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

# print the stoichiometry of first reaction
print rs[0].stoichiometry