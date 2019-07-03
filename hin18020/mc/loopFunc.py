
import subprocess
import os

env = dict(os.environ)
env['LD_LIBRARY_PATH'] = "/Applications/root_v6.14.04/lib"

filelist = "pwglheList.txt"
outputname = "./tmp/output"
otype = ".root"
func_name = "./pythia8Hadronizer"
#func_name = "/Users/tabris/frameLite/hin18020/mc/hadronizer.o"
merge_name = "pwgPythia8.root"

def callfunc(func, arg):
    cmdline = []
    cmdline.append(func)
    cmdline+=arg
    print(cmdline)
    subprocess.call(cmdline, env =env)

def haddmerge():
    cmdline = ["hadd", "-f"]
    cmdline.append(merge_name)
    cmdline.append(outputname+"*"+otype)
    print(cmdline)
    subprocess.call(cmdline)
    #cmdline = ["rm", "-f", outputname+"*"+otype]
    #print(cmdline)
    #subprocess.call(cmdline)

if __name__ == "__main__" :
    istest = 0
    with open(filelist) as f:
       path = f.readlines()
       n = 0
       for p in path:
           if istest and n != 0 : continue
           oname = outputname+"_"+str(n)+otype
           arg = []
           arg.append(p.rstrip())
           arg.append(oname)
       #    if n < 14 : continue
           callfunc(func_name, arg)
           n+=1
    haddmerge()
    """
    """

