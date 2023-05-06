from treeswift import *
from math import log2
AMINOS = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
M62 = [[4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2],
       [0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
       [-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3],
       [-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2],
       [-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3],
       [0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3],
       [-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2],
       [-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1],
       [-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2],
       [-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1],
       [-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1],
       [-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2],
       [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3],
       [-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1],
       [-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2],
       [1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2],
       [0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2],
       [0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1],
       [-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2],
       [-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7]]

BLOSUM62 = dict()
for a in AMINOS:
    BLOSUM62[a] = dict()

N = len(AMINOS)
for i in range(N):
    for j in range(N):
        BLOSUM62[AMINOS[i]][AMINOS[j]] = M62[i][j]

def MSA(sequences, tree, indel_rate):
    '''
    This function computes the MSA for the input sequences
    :param sequences: sequences is a dictionary where keys are species names and values are amino acid sequences
    :param tree: tree is a TreeSwift object
    :param indel_rate: indel_rate is a floating point number giving the relative rate of indel to substitutions
    '''
    aln = dict()

    #function of indel rate
    gap = log2(indel_rate*.67)


    def alignment(ntoalign):
        #  what i L2 is a dicitionary
        align={}
        L1 = aln[ntoalign[0]]
        L2 = aln[ntoalign[1]]
        S = [[None for j in range(len(L2[0])+1)] for i in range(len(L1[0])+1)]
        # fill in axes
        #import pdb; pdb.set_trace()
        # for each letter in first sequence
        S[0][0] = (0,None)
        for i in range(1, len(L1[0])+1): # s1 axis
            #import pdb; pdb.set_trace()
            S[i][0] = (len(L2)*len(L1)*gap*i, 1)
        for j in range(1, len(L2[0])+1): # s2 axis
            S[0][j] = (len(L2)*len(L2)*gap*j, 2)
        
        
        # fill in face
        #import pdb; pdb.set_trace()
        for i in range(1, len(L1[0])+1):
            for j in range(1, len(L2[0])+1):
                # add a gap
                blosum = 0
                for s1 in L1:
                    for s2 in L2:
                        #import pdb; pdb.set_trace()
                        if s1[i-1]  == '-' and s2[j-1] == '-':
                            blosum = blosum + 0
                        elif  s1[i-1]  == '-' or  s2[j-1] == '-':
                            blosum = blosum + gap
                        else:
                            blosum = blosum + BLOSUM62[s1[i-1]][s2[j-1]]
                #import pdb; pdb.set_trace()
                option = [blosum + S[i-1][j-1][0], S[i-1][j][0] + gap*len(L2)*len(L1), S[i][j-1][0] + gap*len(L2)*len(L1)]  
                max_index = max(enumerate(option), key=lambda x: x[1])[0]
                arrow = max_index
                S[i][j] = (max(option), arrow)
        #import pdb; pdb.set_trace()


        # backtracking adapted from given code:
        i = len(L1[0])
        j = len(L2[0])
        alignL1 = ['' for i in range(len(L1))]
        alignL2 = ['' for j in range(len(L2))]

        while arrow is not None:
            if arrow == 1:
                # add gap to L2 strings
                num = 0
                for s1 in L1:  
                    alignL1[num]+=(s1[i-1])
                    num +=1
                num = 0
                for s2 in L2:
                    alignL2[num]+=('-')
                    num += 1
                i-=1
                # add gaps to L1 strings
            elif arrow == 2:
                num = 0
                for s1 in L1:  
                    alignL1[num]+=('-')
                    num +=1
                num = 0
                for s2 in L2:
                    alignL2[num]+=(s2[j-1])
                    num += 1
                j-=1
                
            elif arrow == 0:
                num = 0
                for s1 in L1:  
                    alignL1[num]+=(s1[i-1])
                    num +=1
                num = 0
                for s2 in L2:
                    alignL2[num]+=(s2[j-1])
                    num += 1
                j-=1
                i-=1
            arrow = S[i][j][1]

        #import pdb; pdb.set_trace()
        # invert
        num = 0
        for s1 in alignL1:
            alignL1[num] = s1[::-1]
            num+=1
        num = 0
        for s2 in alignL2:
            alignL2[num] = s2[::-1]
            num+=1
        #import pdb; pdb.set_trace()
        
        #import pdb; pdb.set_trace()

        align
        return alignL1 + alignL2




    for node in tree.traverse_postorder():
        #import pdb; pdb.set_trace()
        if node.is_leaf():
            #import pdb; pdb.set_trace()
            aln[node] = [sequences[node.get_label()]]
            #print(sequences[node.get_label()])
        else:
            #import pdb; pdb.set_trace()
            aln[node] = alignment(node.child_nodes())

    # loop through postorder
    # go from list to dictionary
    
    #import pdb; pdb.set_trace()
    alinlist = aln[tree.root]
    finalalign = {}
    num = 0
    for node in tree.traverse_postorder():            
        if node.is_leaf():
            finalalign[node.get_label()] = alinlist[num]
            num = num +1
        #else: 
         #   import pdb; pdb.set_trace()
          #  finalalign[node.get_label()] = ''.join([finalalign[child.get_label()] for child in node.child_nodes()] )





    #import pdb; pdb.set_trace()

    #raise Exception("NOT IMPLEMENTED ERROR!") # TODO: REPLACE WITH YOUR CODE!
    return finalalign 
