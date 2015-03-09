'''
Created on Mar 9, 2015

@author: mkoenig
'''
#######################################################################
DATA_DIR = "../../../data"
RESULTS_DIR = "../../../results"
VERSION = 2
sql_in = "{}/data.sql".format(DATA_DIR)
sql_out = "{}/data_patched.sql".format(RESULTS_DIR)
#######################################################################

# The tables engine has to be replaced
# ENGINE=MyISAM -> ENGINE=InnoDB

infile = open(sql_in)
outfile = open(sql_out, 'w')

replacements = {'ENGINE=MyISAM':'ENGINE=InnoDB'}

for line in infile:
    for src, target in replacements.iteritems():
        line = line.replace(src, target)
    outfile.write(line)
infile.close()
outfile.close()

print sql_in, '->', sql_out