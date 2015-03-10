'''
Created on Mar 10, 2015

@author: mkoenig
'''
import os
from libsbml import *
from metabolism_settings import RESULTS_DIR, VERSION
from django.core.exceptions import ObjectDoesNotExist

from public.models import Protein, EnzymeParticipant, Entry
from public.models import Reaction as DBReaction 

# Create list of input-output requirements
def calculate_requirements(m):
    from pandas import DataFrame
    print '*** Requirements ***'
    # list
    cols = ['id', 'name', 'readonly', 'readwrite', 'requirement']
    info = []
    
    # species
    for s in m.getListOfSpecies():
        sid = s.getId()
        tokens = sid.split('_')
        # remove the unnecessary prefixes
        wid = "_".join(tokens[1:(len(tokens)-1)])
        cid = tokens[len(tokens)-1] 
        
        try:
            # Check if entry
            e = Entry.objects.get(wid=wid)
            info.append(['_'.join([wid, cid]), e.name, '-', 'x', 'x'])
            
        except ObjectDoesNotExist:
            print 'No wid', sid, wid
        
    # reactions 
    for r in m.getListOfReactions():
        sid = r.getId()
        # remove the unnecessary prefixes
        tokens = sid.split('_')
        # remove the unnecessary prefixes
        wid = "_".join(tokens[1:len(tokens)])
        try:
            # Check if entry
            db_r = DBReaction.objects.get(wid=wid)
            # get the 
            enzyme = db_r.enzyme
            if (enzyme != None):
                protein = enzyme.protein
                # print protein.wid
                info.append([protein.wid, protein.name, 'x', '-', '-'])
            else:
                print 'Reaction has no enzyme', wid, '\t',  db_r.name
            
        except ObjectDoesNotExist:
            print 'No wid', sid, wid
    
    print '#' * 30 
    info_df = DataFrame(info, columns=cols)
    return info_df


if __name__ == '__main__':

    sbml = os.path.join(RESULTS_DIR, "Metabolism_annotated_{}.xml".format(VERSION))
    out_csv = os.path.join(RESULTS_DIR, "requirements.csv".format(VERSION))
    
    doc = readSBML(sbml)
    m = doc.getModel()
    
    info_df = calculate_requirements(m)
    print info_df.head(10)
    info_df.to_csv(out_csv, sep='\t',encoding='UTF8')
    print out_csv


    # for d in info:
    #     print d
     