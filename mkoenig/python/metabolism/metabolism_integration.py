'''
Information for the integation group.
Namely the input-output requirements of the Process Metabolism.

@author: Matthias Koenig
@date: 2015-03-12
'''
import os
from libsbml import *
from metabolism_settings import RESULTS_DIR, VERSION
from django.core.exceptions import ObjectDoesNotExist

from public.models import Protein, EnzymeParticipant, Entry
from public.models import Reaction as DBReaction 

from pandas import DataFrame

# Create list of input-output requirements
def calculate_requirements(m):
    print '*** Requirements ***'
    # list
    cols = ['id', 'model_type', 'name', 'readonly', 'readwrite', 'requirement']
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
            if (cid != 'c'):
                info.append(['__'.join([wid, cid]), e.model_type, e.name, '', 'x', 'x'])
            else: 
                info.append([wid, e.model_type, e.name, '', 'x', 'x'])
            
            
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
                info.append([protein.wid, protein.model_type, protein.name, 'x', '', ''])
            else:
                print 'Reaction has no enzyme', wid, '\t',  db_r.name
            
        except ObjectDoesNotExist:
            print 'No wid', sid, wid
    
    print '#' * 30 
    info_df = DataFrame(info, columns=cols)
    return info_df


if __name__ == '__main__':

    sbml = os.path.join(RESULTS_DIR, "Metabolism_annotated_{}.xml".format(VERSION))
    out_csv = os.path.join(RESULTS_DIR, "Requirements_Metabolism.csv".format(VERSION))
    
    doc = readSBML(sbml)
    m = doc.getModel()
    
    info_df = calculate_requirements(m)
    print info_df.head(10)
    info_df.to_csv(out_csv, sep='\t',encoding='UTF8')
    
    # write csvs
    for name in ['Metabolite', 'ProteinMonomer', 'ProteinComplex']:
        f_csv = os.path.join(RESULTS_DIR, 'Molecules_names_{}_Metabolism.csv'.format(name))
        # write the subsets
        df = info_df[info_df['model_type']==name]
        print df.head()
        df = df[['id', 'name', 'readonly', 'readwrite', 'requirement']]
        df = df.drop_duplicates()
        import csv
        df.to_csv(f_csv, sep=',',encoding='UTF8', index=False, quoting=csv.QUOTE_ALL)
    
        
    # parameters changed  
    f_csv = os.path.join(RESULTS_DIR, 'Parameters_Metabolism.csv'.format(name))
    df = DataFrame()
    df.to_csv(f_csv, sep=',',encoding='UTF8', index=False, quoting=csv.QUOTE_ALL)
