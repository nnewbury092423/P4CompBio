import unittest
from msa.todo import *
from msa.utils import *
from os.path import dirname, realpath
from os import listdir, remove
from subprocess import run
from tempfile import NamedTemporaryFile
from treeswift import *
from func_timeout import func_timeout, FunctionTimedOut 

TIMEOUT = 600 # seconds

test_path = dirname(realpath(__file__))+"/test_data/checking/"
test_cases = ['lowrate','highrate']
test_ids = {'lowrate':['1','2'],'highrate':['1','2']}
rate = {'lowrate':0.01, 'highrate':0.1}
cutoff = {'lowrate':0.95,'highrate':0.85}
threshold = 0.8

class Tests_01(unittest.TestCase):
    def test_01_sanity(self):
        # test reading guidance tree
        tree_file = test_path + "/guidance.nwk"
        try:
            read_tree_newick(tree_file)
        except:
            self.assertTrue(False,msg="Couldn't read the guidance tree! Please check the installation of TreeSwift and make sure the tree file " + tree_file + " exists!")    

        # test reading sample inputs
        for case in test_cases:
            for ID in test_ids[case]:
                test_file = test_path + "/" + case + "/" + ID + ".fas" 
            try:
                with open(test_file,'r') as f:
                    seqs = read_FASTA(f)
            except:        
                    self.assertTrue(False,msg="Couldn't read the sample input file " + test_file + "!")

    def test_02_sanity(self):
        # test running FastSP
        for case in test_cases:
            for ID in test_ids[case]:
                test_file = test_path + "/" + case + "/" + ID + "_TRUE.fas" 
            tempSP = NamedTemporaryFile(delete=False)
            try:
                run(["java","-jar","FastSP.jar","-r",test_file,"-e",test_file,"-o",tempSP.name],stderr=NamedTemporaryFile(mode='w'),check=True)
                with open(tempSP.name,'r') as f:
                    sp = 0
                    md = 0
                    for line in f:
                        if line.startswith("SP-Score"):
                            sp = float(line.strip().split("SP-Score")[-1])
                        if line.startswith("Modeler"):    
                            md = float(line.strip().split("Modeler")[-1])
                    self.assertTrue(sp == 1.0 and md == 1.0,msg="Wrong run of FastSP on " + test_file)
            except:        
                    self.assertTrue(False,msg="Couldn't run FastSP on " + test_file)

            remove(tempSP.name)

    def __run_case_sanity__(self,case,ID):
        sample_in = test_path + "/" + case + "/" + ID + ".fas"
        tree_file = test_path + "/guidance.nwk"
        
        with open(sample_in,'r') as f:
            seqs = read_FASTA(f)
        tree = read_tree_newick(tree_file)    
        try:
            aln = func_timeout(TIMEOUT,MSA,args=(seqs,tree,rate[case],))
        except: 
            self.assertTrue(False,msg="Failed sanity test on " + case + " " + ID + ": Couldn't produce output in time!")
       
        L = None
        for sp in aln:
            L1 = len(aln[sp])
            if L is None:
                L = L1
                self.assertTrue(L == L1, msg="Failed sanity test on " + case + " " + ID + ":Aligned sequences does not have equal length")
            self.assertTrue(aln[sp].replace('-','') == seqs[sp], msg="Failed sanity test on " + case + " " + ID + ": Removing gaps from the aligned sequences does not yield the original sequences")
    
    def __run_case_correctness__(self,case,ID):
        sample_in = test_path + "/" + case + "/" + ID + ".fas"
        sample_ref = test_path + "/" + case + "/" + ID + "_TRUE.fas"
        tree_file = test_path + "/guidance.nwk"        
        
        with open(sample_in,'r') as f:
            seqs = read_FASTA(f)
        tree = read_tree_newick(tree_file)
            
        try:
            aln = func_timeout(TIMEOUT,MSA,args=(seqs,tree,rate[case],))
        except: 
            self.assertTrue(False,msg="Failed correctness test on " + case + " " + ID + ": Couldn't produce output in time!")
        
        tempOut = NamedTemporaryFile(delete=False)
        
        with open(tempOut.name,'w') as f:
            for ID in aln:
                f.write(">%s\n%s\n" % (ID,aln[ID]))
        
        tempSP = NamedTemporaryFile(delete=False)
        run(["java","-jar","FastSP.jar","-r",sample_ref,"-e",tempOut.name,"-o",tempSP.name],stderr=NamedTemporaryFile(mode='w'),check=True)
        
        with open(tempSP.name,'r') as fin:
            for line in fin:
                if line.startswith("SP-Score"):
                    score = float(line.split("SP-Score")[-1].strip())
                    break

        self.assertTrue(score >= threshold, msg = "Failed correctness test on" + case + " " + ID + ": Too low SP-Score. Expect: >=" + str(threshold) + ". Your score: " + str(score))
        if score < cutoff[case]:
            print("Partially passed correctness test on" + case + " " + ID + ": Your SP-Score is " + str(score) + ". Aim at score " + str(cutoff[case]) + " to get full credit") 
        remove(tempOut.name)
        remove(tempSP.name)
    
    def test_03_sanity(self):
        case = 'lowrate'
        ID = '1'
        self.__run_case_sanity__(case,ID)
    
    def test_04_sanity(self):
        case = 'lowrate'
        ID = '2'
        self.__run_case_sanity__(case,ID)
    
    def test_05_sanity(self):
        case = 'highrate'
        ID = '1'
        self.__run_case_sanity__(case,ID)
    
    def test_06_sanity(self):
        case = 'highrate'
        ID = '2'
        self.__run_case_sanity__(case,ID)
   
    def test_07_correctness(self):
        case = 'lowrate'
        ID = '1'
        self.__run_case_correctness__(case,ID)
    
    def test_08_correctness(self):
        case = 'lowrate'
        ID = '2'
        self.__run_case_correctness__(case,ID)

    def test_09_correctness(self):
        case = 'highrate'
        ID = '1'
        self.__run_case_correctness__(case,ID)
    
    def test_10_correctness(self):
        case = 'highrate'
        ID = '2'
        self.__run_case_correctness__(case,ID)
