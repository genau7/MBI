
seqRef1="GAATTC"
seq1="GATTA"

seqRefk="AACACTTTTCAAT"
seqk="ACTTATCA"

seqRef="CCATCGCCATCG"
seq="GCATCGGC"

penalty=-5
k=2
rows=0
cols=0

def getTuplesList():
	tuples=set()
	tuplesDict={}
	for i in range(0, len(seq)-k):
		tuple=seq[i:i+k]
		if tuple not in tuples:
			tuples.add(tuple)
			tuplesDict[tuple]=i
	return tuples,tuplesDict
	
def createDiagonalDict():
	#TODO: verify that
	diagonalNum=rows+cols-2*k+1
	diagonalIndexDict={}
	for row in range(0, rows-k+1):
		diagonalIndexDict[row]=0
	for col in range(1, cols-k+1):
		diagonalIndexDict[-col]=0
	return diagonalIndexDict
	
def calcDiagonalSums():
	cols=len(seqRef)
	rows=len(seq)
	diagonalSum=createDiagonalDict()
	
	scoreMatrix=[[0 for col in range(cols)] for row in range(rows)]

	#TODO boundaries ??(+1)
	for row in range(0, rows-k+1):
		tuple=seq[row:row+k]
		for col in range(0, cols-k+1):
			if(tuple==seqRef[col:col+k]):
				print "incrementing diagonal ", col
				offset=row-col
				diagonalSum[offset]+=1
				scoreMatrix[row][col]=1
				print "  ",tuple, "found at ", row, col
	return diagonalSum, scoreMatrix

def matrixWithTopDiagonals(diagonals, scoreMatrix):
	#get keys for the top ten diagonal sums
	
	#remove entries with score=0
	filteredDiagonals=dict(diagonals)
	for key in diagonals:
		if diagonals[key]==0:
			del filteredDiagonals[key]			
			
	keysForBestDiagonals=sorted(filteredDiagonals, key=filteredDiagonals.__getitem__, reverse=True)
	if len(keysForBestDiagonals)>10:
		keysForBestDiagonals=keysForBestDiagonals[0:11]

	
	#get updated diagonals dictionary
	bestDiagonals=dict(filteredDiagonals)
	print"---then---"
	print bestDiagonals
	for diagKey in filteredDiagonals.keys():
		if diagKey not in keysForBestDiagonals:
			del bestDiagonals[diagKey]
	print "----now---"
	print bestDiagonals
	
	for row in range(0,len(scoreMatrix)):
		for col in range(0, len(scoreMatrix[0])):
			if scoreMatrix[row][col]==1:
				#check if it lays on one of the top 10 diagonals
				diagKey=row-col
				scoreMatrix[row][col]=0 if diagKey not in keysForBestDiagonals else 1
				
	return scoreMatrix, bestDiagonals
	"""
	print "best diags:"
	print keysForBestDiagonals
	print
	for i in range (len(keysForBestDiagonals)):
		print keysForBestDiagonals[i]"""
	
def rescoreRegions(scoreMatrix,blosum,bestDiagonals):

        """ Recore the matrix using blosum matirx and trim notmatching diagonal subregions
        !!!what should be done here?
        - rescore only for m[i][j]!=0 for each good diagonal?
        - rescore for each m[i][j] on each good diagonal?

        """
	diagonalSums=createDiagonalDict()

				
	#iterate over diagonals and score the with blosum matrix
	for diag in bestDiagonals:
		offset=diag
		firstRow=0 if diag<0 else offset
		print "Analizing diag=", diag
		for row in range(firstRow, rows-k-1):
			col=row-offset
			if(col<=cols-k+1):
				#print "row=", row, "col=", col, "offset=", offset, seq[row], seqRef[col]
				if(scoreMatrix[row][col]!=0):
					newValue=blosum[(seq[row], seqRef[col])]
					scoreMatrix[row][col]=newValue
					diag+=newValue
		
	printMatrix(scoreMatrix)
	#find subregions and trim the ones not optimal
	for diag in bestDiagonals:
		regionsSorted=findRegions(diag, scoreMatrix)
		print regionsSorted
		offset=diag
		firstRow=0 if diag<0 else offset
		for row in range(firstRow, rows-k+1):
			if row not in regionsSorted[0][1]:
				col=row-offset
				if(col<=cols-k+1):
					scoreMatrix[row][col]=0
		
	return scoreMatrix,diagonalSums
	
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
			if scoreMatrix[i][j]==1:
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
blosum=readBlosum("blosum.txt")

print "Comparing:"
print "Query=    ", seq, " with"
print "Reference=", seqRef

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


	
