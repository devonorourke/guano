'''
Created on Dec 30, 2015

@author: Toni Westbrok
@contact: anthonyw@wildcats.unh.edu

@summary: Count total occurrence and per-line existence count per query sequence
  1.  If second argument given, query sequences are each line of text file (assumed 8BP)
  2.  If no second argument, query sequences are the first 8BP of each sequence in the fastq file
  3.  Total count indicates the total number of times the query sequence was found
  4.  Line count indicates the number of times the query sequence existed in a read sequence
  5.  Starting tag count indicates the number of times the query sequence was the start of a read sequence
  6.  Only sequences starting with a query (regardless of query source) will be included in results
  
  Usage: devon-count.py <input.fq> [querylist.txt] 

'''

import sys

# Read sequences from fastq
def readSeqs():
    retSeqs = list()
    inFile = open(sys.argv[1], "r")
    
    blockLine = 0
    for currentLine in inFile:
        if blockLine == 1:
            retSeqs.append(currentLine.strip())
            
        blockLine = (blockLine + 1) % 4

    return retSeqs

# Filter seqs for those starting with query values
def filterSeqs(passSeqs, passQueries):
    retSeqs = list()
    
    for currentSeq in passSeqs:
        if len(currentSeq) >= 8:
            if currentSeq[:8] in passQueries:
                retSeqs.append(currentSeq)
                
    return retSeqs

# Import queries from a file
def importQueries():
    retQueries = set()
    inFile = open(sys.argv[2], "r")
    
    for currentLine in inFile:
        retQueries.add(currentLine.strip())
        
    return retQueries

# Parse a unique list of all query strings (first 8 BP)
def getQueries(passSeqs):
    retQueries = set()
    
    for currentSeq in passSeqs:
        if len(currentSeq) >= 8:
            retQueries.add(currentSeq[:8])

    return retQueries

# Report total number of occurrences
def reportTotal(passSeqs, passQueries):
    print("\nTotal Counts:")
    
    for currentQuery in passQueries:
        currentCount = 0
        for currentSeq in passSeqs:
            currentCount += currentSeq.count(currentQuery)
            
        print("{0}: {1}".format(currentQuery, currentCount))

# Report line number of occurrences
def reportLine(passSeqs, passQueries):
    print("\nLine Counts:")
    
    totalCount = 0
    for currentQuery in passQueries:
        currentCount = 0
        for currentSeq in passSeqs:            
            currentCount += currentQuery in currentSeq
            
        print("{0}: {1}".format(currentQuery, currentCount))

# Report starting tag number of occurrences
def reportTag(passSeqs, passQueries):
    print("\nStarting Tag Counts:")
    
    totalCount = 0
    for currentQuery in passQueries:
        currentCount = 0
        for currentSeq in passSeqs:            
            currentCount += currentSeq.startswith(currentQuery)
            
        print("{0}: {1}".format(currentQuery, currentCount))
    
# Get sequences
seqs = readSeqs()

# Read queries from file is specified, otherwise detect from first 8BP
if len(sys.argv[1:]) == 2:
    queries = importQueries()
else:
    queries = getQueries(seqs)
queries = sorted(queries)

# Filter seqs for ones starting with our query values
seqs = filterSeqs(seqs, queries)

# Report counts
reportTotal(seqs, queries)
reportLine(seqs, queries)
reportTag(seqs, queries)
