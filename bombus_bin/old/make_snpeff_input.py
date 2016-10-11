from  __future__ import division
import sys
import os

def task():
	if len(sys.argv) != 3:
		print >> sys.stderr, "Usage Incorrect!\nCorrect Usage: make_snpeff_input.py	InputFile	OutputFile"
		sys.exit()

	inputfile = sys.argv[1]
	output = sys.argv[2]

	rowOutput = "ID\tXP\tLength\n"
	outputFileWriter = open(output, "w")
	outputFileWriter.write(rowOutput)
	outputFileWriter.close()

	id_global = []
	xp_global = []
	length_global = []
	global_count = 0

	with open(inputfile, "r") as inputline:	
		
		for line in inputline:
			
			idd = line.strip().split()[0]		
			xp = line.strip().split()[1] 
			length = line.strip().split()[2]
			
			if (idd not in id_global):
				
				if (global_count == 0):
					global_count = global_count + 1
					
					id_global.append(idd)	
					xp_global.append(xp)
					#print
					length_global.append(length)					

				else : 
					if (len(xp_global) > 1):
						i = length_global.index (max(length_global))
						q = xp_global[i]
					else:
						q = xp_global[0]
					
					rowOutput = id_global[0] + "\t" + max(length_global) + "\t" + q + "\n"
					outputFileWriter = open(output, "a")
					outputFileWriter.write(rowOutput)
					outputFileWriter.close()

					xp_global = []
					id_global = []
					length_global = []
					
					id_global.append(idd)	
					xp_global.append(xp)
					#print
					length_global.append(length)	
					
			elif (idd in id_global):
				xp_global.append(xp)
				length_global.append(length)							


if __name__ == "__main__":
	import sys
	import collections
	from collections import deque
	from collections import defaultdict
	import csv
	import os, os.path
	import time
	import multiprocessing
	import glob	
	task()


