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
# find some reactions
Reaction.objects.filter(name__startswith='A')
# find some genes
Gene.objects.filter(name__contains="2")
