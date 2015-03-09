'''
Created on Mar 9, 2015

@author: mkoenig
'''
import sys
sys.path.append('/home/mkoenig/wholecell-metabolism/mkoenig/python/database/')
sys.path.append('/home/mkoenig/wholecell-metabolism/mkoenig/python/database/wholecellkb')
sys.path.append('/home/mkoenig/wholecell-metabolism/mkoenig/python/database/wholecellkb/wholecellkb')
sys.path.append('/home/mkoenig/wholecell-metabolism/mkoenig/python/database/wholecellkb/public')

# sys.path.append('/home/mkoenig/multiscale-galactose/python')

# has to be overwritten
# os.environ['DJANGO_SETTINGS_MODULE'] = 'wholecellkb.settings'
# import wholcellkb.settings
from public.models import Metabolite, Reaction, Gene


def test():
    # metabolites = Metabolite.objects.all()[:10]
    # print metabolites
    pass


if __name__ == '__main__':
    import sys
    for path in sys.path: print path
    # test()
