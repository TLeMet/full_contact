#!usr/bin/env python

# TLeMet April 3, 2019
# Program for counting contacts between amino acids and lipids
# Please make sure GROMACS is sourced/loaded

import sys
import os
import argparse
from optparse import OptionParser

os.system('mkdir special_contact')

usage = "$ special_contact.py -f *.xtc -s *.gro -l 'list of lipids to be measured'"
parser = OptionParser(usage = usage)

parser.add_option("-f", "--trajectory", action="store", dest="trajectory", default=False)
parser.add_option("-s", "--initial", action="store", dest="initial", default=False)
parser.add_option("-l", "--lipids", action="store", dest="lipids", default=False)

(options, args) = parser.parse_args()

for input_file in [options.trajectory, options.initial, options.lipids]:
	if(input_file==False):
		sys.exit("\nError: Missing input.")

os.system('echo \'q \n\' | gmx make_ndx -f '+options.initial+' -o system.ndx')
os.system('egrep "[\[*]" system.ndx > bottom_text')
groups = list(open('bottom_text','r'))
os.system('rm bottom_text')
groups = [groups[i].strip().replace('[','').replace(']','').strip() for i in range(len(groups))]
index_protein = groups.index('Protein')

options.lipids = [i.strip() for i in options.lipids.upper().split(',')]
lipid_groups = [i for i in range(len(groups)) if groups[i] in options.lipids]

'''
for i in lipid_groups:
	os.system('echo \''+str(index_protein)+' \n '+str(i)+' \n\' | gmx mindist -d 0.8 -f '+options.trajectory+' -s '+options.initial+' -n system.ndx -od special_contact/mindist_Protein_'+str(groups[i])+'.xvg -on full_contact/numcont_Protein_'+str(groups[i])+'.xvg')

membrane = '|'.join(str(lipid_groups).replace('[','').replace(']','').replace(' ','').split(','))
'''

os.system('echo \' splitres 1 \n q \n\' | gmx make_ndx -f '+options.initial+' -o system.ndx -n system.ndx')
os.system('egrep "[\[*]" system.ndx > bottom_text')
groups = list(open('bottom_text','r'))
os.system('rm bottom_text')
groups = [groups[i].strip().replace('[','').replace(']','').strip() for i in range(len(groups))]
#index_membrane = groups.index('Membrane')

aa_groups = [i for i in range(len(groups)) if (('Protein_' in groups[i]) and (not 'Membrane' in groups[i]))]

for i in aa_groups:
	for j in lipid_groups:
		os.system('echo \''+str(j)+' \n '+str(i)+' \n\' | gmx mindist -d 0.8 -f '+options.trajectory+' -s '+options.initial+' -n system.ndx -on special_contact/numcont_'+str(groups[j])+'_'+str(groups[i][8:])+'.xvg')
