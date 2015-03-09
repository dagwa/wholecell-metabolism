'''
Test the object mapping.
@author: mkoenig
@date: 2015-03-09
'''
import sys
from public.models import Metabolite, Reaction, Gene

def test():
    metabolites = Metabolite.objects.all()[:10]
    
    for m in metabolites: 
        print m.pk, m.name, m.charge


if __name__ == '__main__':
    test()
