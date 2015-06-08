
seqRef1="GAATTC"
seq1="GATTA"

seqRefk="AACACTTTTCAAT"
seqk="ACTTATCA"

seqRef="CCATCGCCATCGG"
#seq="GCATCGGC"
seq="ACATCTCCAGCG"



#ktup lenght; should be 4,5 or 6 for DNA, preferably 6
k=2
# some value to filter out diagonals with lowest scores (below cutoff means low score)
cutoff=10 


def getTuplesList(str):
	tuples=set()
	tuplesDict={}
	
	for i in range(0, len(str)-k+1):
		tuple=str[i:i+k]
		
		if tuple not in tuples:
			tuples.add(tuple)
			tuplesDict[tuple]=[]
			
		tuplesDict[tuple].append(i)
	return tuples,tuplesDict
	
def createDiagonalDict():
	#TODO: verify that
	#diagonalNum=rows+cols-2*k+1

	diagonalIndexDict={}
	for row in range(0, rows-k+1):
		diagonalIndexDict[row]=0
	for col in range(1, cols-k+1):
		diagonalIndexDict[-col]=0
	return diagonalIndexDict
	
def createDiagonalDictFrom(bestDiagonals):
	newDict=dict(bestDiagonals)
	for diag in bestDiagonals:
		newDict[diag]=0 
	return newDict	
	
def createHotspotsDict():
	#TODO: verify that
	#diagonalNum=rows+cols-2*k+1

	diagonalIndexDict={}
	for row in range(0, rows-k+1):
		diagonalIndexDict[row]=[]
	for col in range(1, cols-k+1):
		diagonalIndexDict[-col]=[]
	return diagonalIndexDict
	
def calcDiagonalSums(tuplesRef, tuplesRefDict):
	diagonalSums=createDiagonalDict()
	hotspotRows=createHotspotsDict()
	#scoreMatrix=[[0 for col in range(cols)] for row in range(rows)]

	for row in range(0, rows-k+1):
		tuple=seq[row:row+k]	
		if tuple in tuplesRef:
			
			for col in tuplesRefDict[tuple]:
				offset=row-col
				diagonalSums[offset]+=1
				hotspotRows[offset].append(row)
	
	#get rid of entries with sum==0
	answer=dict(diagonalSums)
	for diag in diagonalSums:
		if diagonalSums[diag]==0:
			del answer[diag]
			del hotspotRows[diag]
	
	print"-------------------"
	print "diagonalSums="
	print answer	
	print "-------------------"
	print "hotspotRows"
	print hotspotRows
	print"-----------------------"
	print "tuplesRefDict (hotspotCols):"
	print tuplesRefDict
	return answer, hotspotRows

def scoreDiagonals(diagonalSumDict, hotspotRows):
	"""
	In order to evaluate each diagonal run, FASTA gives each hot spot a positive score,
	and the space between consecutive hot spots in a run is given a negative score that
	decreases with the increasing distance. The score of the diagonal run is the sum of the
	hot spots scores and the interspot scores. FASTA finds the 10 highest scoring diagonal
	runs under this evaluating scheme.

	"""
	rescoredDiagonals=dict(diagonalSumDict)
	reward=20 #this should be some kind of positive value

	for diag in diagonalSumDict:
		print "diag=",diag
		
		firstRow=0 if diag<=0 else diag	
		gapPenalty=-reward/2 #this should be some kind of negative value
		interspotPenaltySum=0
		sum=0
		
		#for each row, score the cell and update diagonalSums
		for row in range(firstRow,rows-k+1):
			col=row-diag
			
			# check if col within boundaries
			if(col<cols):
				hotspot=row in hotspotRows[diag]
				
				#if for this col there is a gap between two hotspots: 
				#(sum>0 indicates at least 1 hotspot was found before)
				if(hotspot==False and sum>0):
					#sum penalties for each consecutive gap
					interspotPenaltySum+=gapPenalty
					
					#make sure penalty stays a negativ value (in case there were too many consecutive gaps)
					if(gapPenalty<0):
						#each consecutive gap has less of a penalty
						gapPenalty+=1
										
				elif hotspot:
					#sum the reward and add the penalty for the gaps between this hotspot and the previous one
					sum=sum+reward+interspotPenaltySum
					
					#since a hotspot was found, reset the interspotPenaltySum and revert to the initial gapPenalty value
					interspotPenaltySum=0
					gapPenalty=-reward/2
				#print "row=", row, "sum=", sum, "gapPenalty=", interspotPenaltySum
		
		rescoredDiagonals[diag]=sum
		sum=0
	
	if(len(rescoredDiagonals)>10):
		rescoredDiagonals=getTopDiagonals(getTopDiagonals)
	return rescoredDiagonals

def getTopDiagonals(diagonals):	
	#get keys for the top ten diagonal sums			
	keysForBestDiagonals=sorted(diagonals, key=diagonals.__getitem__, reverse=True)
	keysForBestDiagonals=keysForBestDiagonals[0:11]

	#get updated diagonals dictionary
	bestDiagonals=dict(diagonals)
	print"---then---"
	print bestDiagonals
	for diagKey in diagonals.keys():
		if diagKey not in keysForBestDiagonals:
			del bestDiagonals[diagKey]
	
	print "----now---"
	print "best diagonals:\n", bestDiagonals
				
	return bestDiagonals
	"""
	print "best diags:"
	print keysForBestDiagonals
	print
	for i in range (len(keysForBestDiagonals)):
		print keysForBestDiagonals[i]"""	


def rescoreDiagonals(blosum,bestDiagonals,hotspotRows):

	rescoredDiagonals=createDiagonalDictFrom(bestDiagonals)
	print rescoredDiagonals			
	#iterate over diagonals and score the with blosum matrix
	for diag in bestDiagonals:
		for row in hotspotRows:
			col=row-diag
			rescoredDiagonals[diag]+=blosum[(seq[row], seqRef[col])]
		
	
	#remove diagonals with scores below a cutoff threshold
	bestRescoredDiagonals=dict(rescoredDiagonals)
	
	print rescoredDiagonals
	
	for diag in rescoredDiagonals:
		if rescoredDiagonals[diag]<cutoff:
			del bestRescoredDiagonals[diag]
			
	print
	print bestRescoredDiagonals
	return bestRescoredDiagonals

	
def createMatrixForDots(diagonalSumDict, hotspotRows):
	"""
	just for debugging. drawing dot matrix is unnecessary since it takes too much memory to remember all values
	"""
	scoreMatrix=[[0 for col in range(cols)] for row in range(rows)]
	for diag in diagonalSumDict:
		firstRow=0 if diag<=0 else diag
		for row in range(firstRow,rows-k+1):
			col=row-diag
			if(col<cols):
				#print "row=", row, "col=", col
				if(row in hotspotRows[diag]):
					scoreMatrix[row][col]=1
				else:
					scoreMatrix[row][col]=0
	return scoreMatrix

	
def findRegions(diag, scoreMatrix):
	regionRows=[]
	regions=[]
	offset=diag
	firstRow=0 if diag<0 else offset
	print "Analizing diag=", diag
	sum=0
	print "rang=", rows-k+1
	for row in range(firstRow, rows-k+1):
		col=row-offset
		#print "rangeCol=",  col
		if(col<=cols-k+1):
			print "row=",row, "offset=",offset,
			if(scoreMatrix[row][col]!=0):
				sum+=scoreMatrix[row][col]
				regionRows.append(row)
				if(row==rows-k and sum>0):
					regions.append((sum, regionRows[:]))
			else:
				if(sum>0):
					regions.append((sum, regionRows[:]))
				sum=0
				regionRows[:]=[]
			print "col=",col, "val=",scoreMatrix[row][col],"sum=",sum
		

	return sorted(regions, key=lambda tup: tup[0], reverse=True)	
		
# Create an empty matrix
def create_matrix(m, n):
    return [[0]*n for _ in xrange(m)]

# Read the BLOSUM50 table
def readBlosum(fname):

    #dictionary holding pairs of nucleotides and their BLOSUM matrix value, e.g.: ('T', 'A'): -4
    d = {}
    lines = open(fname, "rt").readlines()
    alpha = lines[0].rstrip('\n\r').split()
    assert(len(alpha) == len(lines)-1)
    for r in lines[1:]:
        r = r.rstrip('\n\r').split()
        a1 = r[0]
        for a2, score in zip(alpha, r[1:]):
            d[(a1, a2)] = int(score)
    return d

def needlemanWunsch(seqVertical, seqHorizontal, blosum, penalty):
    rows = len(seqVertical)+1 
    cols = len(seqHorizontal)+1 
    
    F = create_matrix(rows, cols)
 
    for i in range(0, rows):
        F[i][0] = i * penalty
    for j in range(0, cols):
        F[0][j] = j * penalty
 
    for i in range(1, rows):
        for j in range(1, cols):
            match = F[i-1][j-1] + blosum[(seqVertical[i-1], seqHorizontal[j-1])]
            delete = F[i-1][j] + penalty
            insert = F[i][j-1] + penalty
            F[i][j] = max(match, delete, insert)
 
    return F


def traceback(scoreMatrix, startPos, blosum):
    '''Find the optimal path through the matrix representing the alignment.

    Starting from the best position (bottom right of a path), trace the whole path
    back up (top-left corner), finding thus the best local alignment. Each step of the path (matrix cell)
    corresponds to either a gap in a sequence (or both sequences) or a match/mismatch in the following way:
        diagonal (i-1, j-1) - match/mismatch
        up       (i-1, j  ) - gap in sequence 1
        left     (i  , j-1) - gap in sequence 2

    A step that should be taken is the one that leads to the predecessor cell
    '''

    END, DIAG, UP, LEFT = range(4)
    alignedSeq = []
    alignedSeqRef = []
    i, j         = startPos
    step         = nextStep(scoreMatrix, i, j,blosum)
    
    while step != END:
        if step == DIAG:
            alignedSeq.append(seq[i - 1])
            alignedSeqRef.append(seqRef[j - 1])
            i -= 1
            j -= 1
        elif step == UP:
            alignedSeq.append(seq[i - 1])
            alignedSeqRef.append('-')
            i -= 1
        else:
            alignedSeq.append('-')
            alignedSeqRef.append(seqRef[j - 1])
            j -= 1
       
        step = nextStep(scoreMatrix, i, j,blosum)
       
    return ''.join(reversed(alignedSeq)), ''.join(reversed(alignedSeqRef))


def nextStep(scoreMatrix, i, j, blosum):
    score=scoreMatrix[i][j]
    diag = scoreMatrix[i - 1][j - 1]
    up   = scoreMatrix[i - 1][j]
    left = scoreMatrix[i][j - 1]

    similarity=blosum[(seq[i-1], seqRef[j-1])]

    if(score==diag+similarity):
        return 1
    if (score==up+penalty):
        return 2
    if (score==left+penalty):
        return 3

    return 0

def createAlignmentString(alignedSeq, alignedSeqRef):
    '''Construct a special string showing identities, gaps, and mismatches.

    This string is printed between the two aligned sequences and shows the
    identities (|), gaps (-), and mismatches (:). As the string is constructed,
    it also counts number of identities, gaps, and mismatches and returns the
    counts along with the alignment string.

    AAGGATGCCTCAAATCGATCT-TTTTCTTGG-
    ::||::::::||:|::::::: |:  :||:|   <-- alignment string
    CTGGTACTTGCAGAGAAGGGGGTA--ATTTGG
    '''
    # Build the string as a list of characters to avoid costly string
    # concatenation.
    idents, gaps, mismatches = 0, 0, 0
    alignmentString = []
    
    for base1, base2 in zip(alignedSeq, alignedSeqRef):
        if base1 == base2:
            alignmentString.append('|')
            idents += 1
        elif '-' in (base1, base2):
            alignmentString.append(' ')
            gaps += 1
        else:
            alignmentString.append(':')
            mismatches += 1

    return ''.join(alignmentString), idents, gaps, mismatches


def print_matrix(matrix):
    '''Print the scoring matrix.

    ex:
    0   0   0   0   0   0
    0   2   1   2   1   2
    0   1   1   1   1   1
    0   0   3   2   3   2
    0   2   2   5   4   5
    0   1   4   4   7   6
    '''
    for row in matrix:
        for col in row:
            print('{0:>4}'.format(col)),
        print

		
def printDotMatrix(matrix):
	#print ref sequence's nukleotides
	print
	print" ",
	for n in range(cols):
		print seqRef[n], 
	print
	#print sequence's nukleotides and "o" if a dot should be printed
	for i in range(0,rows):
		print seq[i],
		for j in range(0,cols):
			if matrix[i][j]==1:
				print "o",
			else:
				print " ",
		print

def printMatrix(matrix):
	#print ref sequence's nukleotides
	print
	print"    ",
	for n in range(cols):
		print seqRef[n], "  ",
	print
	#print sequence's nukleotides and "o" if a dot should be printed
	i=0
	for row in matrix:
		print seq[i],
		i+=1
		for col in row:
			if col!=0:
				print ('{0:>4}'.format(col)),
			else:
				print "    ",
		print

		
print "---------FASTA-----\n"
rows=len(seq)
cols=len(seqRef)
print "Comparing:"
print "Query=    ", seq, " with"
print "Reference=", seqRef
blosum=readBlosum("blosum.txt")#

# 1.identify common k-words between I (seq) and J(seqRef)
print "\n=======================================\n                STEP1\n=======================================\n"
tuplesRef,tuplesRefDict=getTuplesList(seqRef)
diagonalSumsDict, hotspotRows=calcDiagonalSums(tuplesRef,tuplesRefDict)


print 
print diagonalSumsDict

matrix=createMatrixForDots(diagonalSumsDict, hotspotRows)
printDotMatrix(matrix)

# 2a. Score diagonals with k-word matches and identify 10 best diagonals
print "\n=======================================\n                STEP2\n=======================================\n"
betsTenDiagonals=scoreDiagonals(diagonalSumsDict, hotspotRows)
print "bestTenDiagonals:\n", betsTenDiagonals


# 3. Rescore initial regions with a substitution score matrix and get best 10 subregions
print "\n=======================================\n                STEP3\n=======================================\n"
rescoredDiagonals=rescoreDiagonals(blosum, betsTenDiagonals, hotspotRows)
print "rescoredDiagonals:\n", rescoredDiagonals











""""




#matrix=needlemanWunsch(seq,seqRef, d, penalty)
#print_matrix(matrix)
#rightBottomCell=(len(seq), len(seqRef))
#seqAligned, seqRefAligned = traceback(matrix, rightBottomCell, d)
#print seqRefAligned
#print seqAligned

# 1.identify common k-words between I (seq) and J(seqRef)
# 2a. Score diagonals with k-word matches
diagonalSums, scoreMatrix=calcDiagonalSums()

#print dot mattrix
printDotMatrix(scoreMatrix)
	
print
printMatrix(scoreMatrix)
	
print "diagonal sums:\n"
print diagonalSums

	
# 2b. identify 10 best diagonals	
filteredScoreMatrix,bestDiagonals=matrixWithTopDiagonals(diagonalSums, scoreMatrix)


# 3. Rescore initial regions with a substitution score matrix and get best 10 subregions
scoreMatrix, diagonalSums =rescoreRegions(scoreMatrix,blosum, bestDiagonals)
#Non-matching ends of the diagonal are trimmed?????
#filteredScoreMatrix=matrixWithTopDiagonals(diagonalSums, scoreMatrix)
printMatrix(scoreMatrix)
print


# 4. Join initial regions using gaps, penalise for gaps

# 5. Perform dynamic programming to find final alignments

	"""
	
	
	
	
	
"""
	The stages in the FASTA algorithm are as follows:
1. We specify an integer parameter called ktup (short for k respective tuples), and we look
for ktup-length matching substrings of the two strings. The standard recommended
ktup values are six for DNA sequence matching and two for protein sequence matching.
The matching ktup-length substrings are referred to as hot spots. Consecutive hot spots
are located along the dynamic programming matrix diagonals. This stage can be done
effciently by using a lookup table or a hash to store all the ktup-length substrings from
one string, and then search the table with the ktup-length substrings from the other
string.

2. In this stage we wish to find the 10 best diagonal runs of hot spots in the matrix. A
diagonal run is a sequence of nearby hot spots on the same diagonal (not necessarily
adjacent along the diagonal, i.e., spaces between these hot spots are allowed). A run
need not contain all the hot spots on its diagonal, and a diagonal may contain more
than one of the 10 best runs we find.
In order to evaluate the diagonal runs, FASTA gives each hot spot a positive score,
and the space between consecutive hot spots in a run is given a negative score that
decreases with the increasing distance. The score of the diagonal run is the sum of the
hot spots scores and the interspot scores. FASTA finds the 10 highest scoring diagonal
runs under this evaluating scheme.

3. A diagonal run species a pair of aligned substrings. The alignment is composed of
matches (the hot spots) and mismatches (from the interspot regions), but it does not
contain any indels because it is derived from a single diagonal. We next evaluate
the runs using an amino acid (or nucleotide) substitution matrix, and pick the best
scoring run. The single best subalignment found in this stage is called init1. Apart
from computing init1, a ltration is performed and we discard of the diagonal runs
achieving relatively low scores.

4. Until now we essentially did not allow any indels in the subalignments. We now try to
combine \good" diagonal runs from close diagonals, thus achieving a subalignment with
indels allowed. We take \good" subalignments from the previous stage (subalignments
whose score is above some specied cuto) and attempt to combine them into a single
larger high-scoring alignment that allows some spaces. This can be done in the following
way:
We construct a directed weighted graph whose vertices are the subalignments found
in the previous stage, and the weight in each vertex is the score found in the previous
stage of the subalignment it represents. Next, we extend an edge from vertex u to
vertex v if the subalignment represented by v starts at a lower row than where the
3.3. BLAST - BASIC LOCAL ALIGNMENT SEARCH TOOL 3
subalignment represented by v ends. We give the edge a negative weight which depends
on the number of gaps that would be created by aligning according to subalignment v
followed by subalignment u. Essentially, FASTA then nds a maximum weight path in
this graph. The selected alignment species a single local alignment between the two
strings. The best alignment found in this stage is marked initn. As in the previous
stage, we discard alignments with relatively low score.

5. In this step FASTA computes an alternative local alignment score, in addition to initn. Recall that init1 denes a diagonal segment in the dynamic programming matrix. We
consider a narrow diagonal band in the matrix, centered along this segment. We observe
that it is highly likely that the best alignment path between the init1 substrings, lies
within the subtable dened by the band. We assume this is the case and compute
the optimal local alignment in this band, using the ordinary dynamic programming
algorithm. Assuming that the best local alignment is indeed within the dened band,
the local alignment algorithm essentially merges diagonal runs found in the previous
stages to achieve a local alignment which may contain indels. The band width is
dependent on the ktup choice. The best local alignment computed in this stage is
called opt.

6. In the last stage, the database sequences are ranked according to initn scores or opt
scores, and the full dynamic programming algorithm is used to align the query sequence
against each of the highest ranking result sequences.

"""
