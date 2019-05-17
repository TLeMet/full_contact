#!usr/bin/env python

#TLeMet, May 17th 2019
#Script for plotting special_contact.py results

import sys
import os
import argparse
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def CurveSmoothing(curve, degree=1):
    smoothed = curve*0
    for i in xrange(0,len(smoothed),1):
        if (i<degree):
            smoothed[i] = np.mean(curve[0:i+degree+1])
        elif (i>len(smoothed)-degree-1):
            smoothed[i] = np.mean(curve[i-degree:len(smoothed)])
        else:
            smoothed[i] = np.mean(curve[i-degree:i+degree+1])
    return smoothed

usage = "$ full_contact_plot.py -l 'lipid,types'"
parser = OptionParser(usage = usage)

parser.add_option("-l", "--lipids", action="store", dest="lipids", default=False)

(options, args) = parser.parse_args()

if(options.types==False):
		sys.exit("\nError: Missing input.")

os.system('ls > bottom_text')
files = list(open('bottom_text','r'))
os.system('rm bottom_text')
files = [files[i].replace('\n','') for i in range(len(files))]

options.lipids = [i.strip() for i in options.lipids.upper().split(',')]
lipid_types = [i for i in range(len(groups)) if groups[i] in options.lipids]

data = {}
for f in files:
	for l in lipid_types:
			if ('numcont' in f) and (l in f):
			data[f] = np.transpose(np.genfromtxt([r for r in open(f,'r').readlines() if not r[0] in ('#', '@')]))

for i in data.keys():
	data[i][0] = [float(x)/1000 for x in data[i][0]]

data['prot_z'] = np.transpose(np.genfromtxt([r for r in open('prot_z.xvg','r').readlines() if not r[0] in ('#', '@')]))
data['membrane_z'] = np.transpose(np.genfromtxt([r for r in open('membrane_z.xvg','r').readlines() if not r[0] in ('#', '@')]))
data['dist'] = [abs(data['prot_z'][i]-data['membrane_z'][i]) for i in range(len(data['prot_z']))]


aa = [x[8:-6] for x in data.keys() if ('numcont' in x) and ('Protein' in x)]
print(aa)
'''
for i in lipid_types:


	fig,ax1=plt.subplots()
	
	for aa 
		ax1.plot(data[i][0], CurveSmoothing(data[i][1],10), marker=',', label=i[16:-4]+', avg='+'{:5.4f}'.format(np.mean(data[i][1])))
	ax2 = ax1.twinx()
	ax2.plot(data['prot_z'][0]/1000, CurveSmoothing(data['dist'][1],100), marker=',', linestyle='dotted', label='distance')

	plt.title('Protein contact with each lipid type')
	ax1.set_xlabel('time (ns)')
	ax1.set_ylabel('average number of contacts per lipid molecule')
	ax2.set_ylabel('distance between protein and membrane COMs (nm)')
	handles1, labels1 = ax1.get_legend_handles_labels()
	handles2, labels2 = ax2.get_legend_handles_labels()
	handles, labels = handles1+handles2, labels1+labels2
	fig.legend(handles, labels, loc=2, bbox_to_anchor=(1.1,1), borderaxespad=0)

	plt.savefig('Specific_lipid_Membrane_contacts_Z.png', format='png', bbox_inches='tight', pad_inches=.5, dpi=300)
	plt.clf()

	'''