# a script for aligning models over a 
# certain region in chimerax
# open chimera x, load all sturctures
# then go to tools -> general -> shell

from chimerax.core.commands import run
for s in range(2, 102):
# change 102 above to however many models you have
	string = 'align #'+str(s)':11-37 #1:11-37
	# align all others to the first one based 
	# on overlap in spatial coordinates from 
	# residues 11-37
	run(session, string)

