
import os
import subprocess

#listf = "list_5TeV2015pp_mixingTree.txt"
#listf = "list_bjetPythia8_mixingTree.txt"
#listf = "list_PYTHIA8_bjet_mixingTree.txt"
#listf = "list_Herwig_dijet100to200_mixingTree.txt"
listf = "list_PYTHIA8_dijet_mixingTree.txt"
startline = 0

execf = 'run_Pythia8_dijet.C'

listfs = open(listf).readlines()
counter = 0
outfstr0 = './tmp/output'
for fs in listfs:	
	counter+=1
	if counter < startline+1: continue
	outfstr = outfstr0+'_'+str(counter)+'.root'
        cmdline = ['root', '-l', '-b', '-q', execf+'(\"'+fs.rstrip()+'\", \"'+outfstr+'\")']
	print str(cmdline[4])
	subprocess.call(cmdline)

