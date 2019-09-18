
import os
import sys
import subprocess
"""
using to merge batch files in lpc eos
usage: python autoMerge path_to_target_files
"""
workdir = "/home/wang4187/utility/"
if workdir[-1] != "/": workdir= workdir+"/"

eosdir_prefix = "/mnt/hadoop/store/user/wangx/"

def merge(path):
	outputs = subprocess.check_output(["ls","-u", path])
	plist = outputs.split("\n")
	doit = raw_input("do you want to merge all .root files ? [y: continue]")
	if doit != "y" : return;
	if not os.path.exists(workdir+"tmp/"):
		os.makedirs(workdir+"tmp/")
	ntotal = 0
	flist = []
	if path[-1] != "/": path = path+"/"
	for line in plist:
		if line[-5: ] == ".root" :
			ntotal = ntotal+1
			flist.append(path+line)
	print("total files: ",ntotal)
	print(flist)
	fname = raw_input("please add output file name:")
	commandstring = ["hadd", "-f", workdir+"tmp/"+fname] + flist
	subprocess.call(commandstring)
	subprocess.call(["mv", workdir+"tmp/"+fname, workdir])
#	subprocess.call(["rm", workdir+"tmp/*"])
	return
	


start = False
if len(sys.argv) > 1 : 
	datadir = sys.argv[1]
	path = eosdir_prefix+datadir
	start = True
else :  path = eosdir_prefix

subprocess.call(["ls", path])


if start : merge(path);
