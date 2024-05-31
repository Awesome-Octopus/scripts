#!/usr/bin/python
# -*- coding: utf-8 -*-
from argparse import ArgumentParser

""" 
Hey ! If you're here, you must be interested in what's inside this script.
Basically, it removes all the lines starting with # (which are the include statements)
and puts a pre-defined stack of includes instead.

The whole ArgumentParser stuff is what handles the -f and -o from the command line, 
for inputs and outputs.
"""

parser = ArgumentParser(description=""" Clean the includes in the topology
Used within the INSANE tutorial""")

# Named arguments
parser.add_argument("-f", "--file", help="The name of the input topology file", required=True)
parser.add_argument("-o", "--output", help="The name of the output topology file", required=True)
parser.add_argument("-i", "--insert", help="The names of molecules from itp files", required=True)
args = parser.parse_args()

# Input / Output
topolname = args.file
out_topolname = args.output

######################## The script ########################

# It's straightforward: load the topology (as list of lines)
topology = open(topolname, "r").readlines()
insert = open(args.insert, "r").readlines()

# Remove the lines that start with #
#topology = [x for x in topology if x.startswith("#") == False]
ndx = None
for i in range(0, len(topology)):
	if 'Protein ' in topology[i] and 'INSANE' not in topology[i]:
		ndx = i
	
if ndx is not None:
	a = topology[0:ndx]
	b = topology[ndx+1:]
	for i in insert:
		a.append(i)
	for i in b:
		a.append(i)
else:
	a = topology


topology = "".join(a)



# Open the output
out_ = open(out_topolname, "w")

# Just write everything into it
out_.write(topology)

# Close the output
out_.close()






