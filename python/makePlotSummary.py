
import sys
from os import listdir
from os.path import isfile, join
import subprocess

class texFile :
    def __init__(self):
        self.fstr = ""
        self.fname = ""
    def article_file_header(self):
        self.fstr = "\\documentclass[a4paper]{article}\n\\usepackage{graphicx}\n\\usepackage{listings}\n\\usepackage{xcolor}\n\\usepackage{enumerate}\n\\usepackage{indentfirst}\n\\usepackage{fancyhdr}\n\\usepackage{hyphenat}\n\\usepackage{eufrak} % Gothic font\n\\usepackage{amsmath,amssymb,dsfont}\n\\usepackage{multicol}\n\\usepackage{balance}\n\\usepackage{subfigure}\n\\usepackage{booktabs}\n\\usepackage{algorithm}\n\\usepackage{algpseudocode}% an improvement from algorithmicx for algorithmic\n\\usepackage[top=2.54cm,bottom=2.54cm,left=2.5cm,right=2.5cm]{geometry} % a4paper standard\n\\begin{document}\n"

    def end_file_string(self):
        self.fstr += "\n\\end{document} "
    def print_string(self):
        print(self.fstr)
    def dump_file(self, name):
        self.end_file_string()
        self.fname = name
        with open(name, 'w') as wf : 
            wf.write(self.fstr)
    def add_plot_string(self, caption, path):
        self.fstr+="\n\\begin{figure}\n\\begin{center}\n\\includegraphics[width=1\\textwidth]{"+path+"}\n"
        self.fstr+="\\caption{"+caption+"}\n"
        self.fstr+="\\end{center}\n\\end{figure}\n\n"

    def tex_compile(self, inputf = ""):
        if inputf == "" : inputf = self.fname
        subprocess.call(['pdflatex', inputf])

    def remove_letter(self, string):
        string.replace("_", " ")
        return string

    def add_eps2pdf_dir(self):
        self.fstr = "\\epstopdfsetup{outdir=./}"

    def make_plot_summary(self, name, form, dataset, mypath):
        self.article_file_header()
        if "eps" in form : add_eps2pdf_dir()
        pics = sorted(listdir(mypath))
        files = [f for f in pics if isfile(join(mypath, f))]
        clearcounter = 0
        for f in files:
            if form not in f : continue
            if f.split("_")[0] != dataset : continue
            caption = f.replace("_", " ")
            self.add_plot_string(caption, mypath+f)
            clearcounter += 1
            if clearcounter == 10 :
                clearcounter = 0
                self.fstr += "\\clearpage\n"
        self.dump_file(name+"_"+dataset+".tex")
        print("tex file generated for set: "+dataset+": "+name+"_"+dataset+".tex")



if __name__ == "__main__" :
    ff = texFile()
    #ff.make_plot_summary("../output/doc/summary_step2", ".pdf", "pp5TeVjet80", "/Users/tabris/frameLite/output/step2/")
    #ff.make_plot_summary("../output/doc/summary_step2", ".pdf", "bjetMC", "/Users/tabris/frameLite/output/step2/")
    ff.make_plot_summary("../output/doc/summary_step2", ".pdf", "dijetMC", "/Users/tabris/frameLite/output/step2/")
    #ff.make_plot_summary("../output/doc/summary_step3", ".pdf", "bjtc", "/Users/tabris/frameLite/output/step3/")
    #ff.make_plot_summary("../output/doc/summary_esm_step2", ".pdf", "pp5TeVjet80", "/Users/tabris/frameLite/output/Esc_step2/")
    #ff.make_plot_summary("summary_step2", ".pdf", "bjetMC", "/Users/tabris/frameLite/output/step2/")
    #ff.article_file_header()
    #ff.add_plot_string("inclJetReco effect rawSig ratio", "inclJetReco_effect_rawSig_ratio-eps-converted-to.pdf")
    #print(ff.remove_letter("inclJetReco_effect_rawSig_ratio-eps-converted-to"))
    #ff.dump_file("test.tex")
    #ff.tex_compile()
    #ff.print_string()
