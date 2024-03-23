#!/usr/bin/env python3
#
# Quintet Compressor
# Osteoclave
# 2012-02-04
#
# The compression format is described in the decompressor.
#
# Modified by Tranquilite
#
# Now creates ideal compression chains by utilizing graph theory.
# Marginally beats the greedy longest match compression.
# Also use KMP search insted of brute force to keep compression times low.
# Also switched from bitstring to bitarray since it was a bit faster.
#
# This code uses python-bitarray and networkx:
# https://pypi.org/project/bitarray/
# https://pypi.org/project/networkx/

import os
import sys
from bitarray import bitarray, util
import networkx as nx

# Define some useful constants.
SEARCH_LOG2 = 8
SEARCH_SIZE = 2 ** SEARCH_LOG2
LOOKAHEAD_LOG2 = 4
LOOKAHEAD_SIZE = 2 ** LOOKAHEAD_LOG2
BIT_PASTCOPY = 0
BIT_LITERAL = 1
WEIGHT_PASTCOPY = 13
WEIGHT_LITERAL = 9
MAX_UNCODED = 1
MAX_CODED = LOOKAHEAD_SIZE + MAX_UNCODED
MIN_CODED = 2

def compress(inBytes):
    # Prepare the memory buffer.
    inBuffer = bytearray(SEARCH_SIZE + len(inBytes))
    inBuffer[:SEARCH_SIZE] = [0x20] * SEARCH_SIZE
    inBuffer[SEARCH_SIZE:] = inBytes

    # Create Graph for storing candidate compresion chain
    dag = nx.DiGraph()

    # Main compression loop.
    for currentIndex in range(SEARCH_SIZE, SEARCH_SIZE + len(inBytes)):
        inputPos = currentIndex - SEARCH_SIZE
        bestIndex = 0
        bestLength = 0

        #Maybe Binary search would be better?
        bestIndex, bestLength = KMPSearch(inBuffer, currentIndex-SEARCH_SIZE, currentIndex )

        # Create vertices in our Dag for every possible valid encoding from our current position to a future position
        # and give it a weight which corresponds to the length of the encoding in bits (9 or 13)

        # For some reason, the decompressor expects the pastcopy
        # source values to be offset by 0xEF. I have no idea why.
        bestIndex = (bestIndex + 0xEF) & 0xFF

        # Create verices for matches of all valid match lengths (if any)
        for length in range(MIN_CODED, bestLength+1):
            #bits = bitstring.Bits(f'uint1={BIT_PASTCOPY}, uint{SEARCH_LOG2}={bestIndex}, uint{LOOKAHEAD_LOG2}={length - 2}')
            bits = bitarray([BIT_PASTCOPY]) + util.int2ba(bestIndex, SEARCH_LOG2) + util.int2ba(length - MIN_CODED, LOOKAHEAD_LOG2)
            dag.add_edge(inputPos, inputPos+length, weight=WEIGHT_PASTCOPY, bitstring=bits)

        # Always add the Literal case as an edge to the next node.
        #bits = bitstring.Bits(f'uint1={BIT_LITERAL}, uint8={inBuffer[currentIndex]}')
        bits = bitarray([BIT_LITERAL]) + util.int2ba(inBuffer[currentIndex], 8)
        dag.add_edge(inputPos, inputPos+1, weight=WEIGHT_LITERAL, bitstring=bits)

    # Prepare for compression.
    #output = bitstring.BitArray(uintle=len(inBytes), length=16)
    outputLength = util.int2ba(len(inBytes), 16, 'little')
    output = bitarray()

    # Compute shortest path and compile output
    nodes = nx.shortest_paths.dijkstra_path(dag, 0, len(inBytes))
    edges = nx.utils.pairwise(nodes)

    for u,v in edges:
        #output.append(dag[u][v]["bitstring"])
        output += dag[u][v]["bitstring"]

    # Return the compressed data.
    return outputLength.tobytes() + output.tobytes()

def BruteForceSearch(inBuffer, currentIndex):
    bestLength = 0
    bestIndex = 0
    # Don't compare past the end of the lookahead buffer.
    # Don't compare past the end of the memory buffer.
    compareLimit = min(
        MAX_CODED,
        len(inBuffer) - currentIndex
    )
    
    # Look for a match in the search buffer. (Brute force)
    for i in range(SEARCH_SIZE):
        # Compare the search buffer to the lookahead buffer.
        # Count how many sequential bytes match (possibly zero).
        currentLength = 0
        for j in range(compareLimit):
            if inBuffer[currentIndex - SEARCH_SIZE + i + j] == inBuffer[currentIndex + j]:
                currentLength += 1
            else:
                break
    
        # Keep track of the largest match we've seen.
        if currentLength > bestLength:
            bestIndex = currentIndex - SEARCH_SIZE + i
            bestLength = currentLength

    return bestIndex, bestLength



# KMP Algorithm
# https://en.wikipedia.org/wiki/Knuth%E2%80%93Morris%E2%80%93Pratt_algorithm
def KMPSearch(buffer, windowStart: int, uncodedStart: int):

    # The uncoded lookahead which we want to find the longest prefix match in our lookbehind ring buffer
    uncoded = buffer[uncodedStart : uncodedStart + MAX_CODED] # I am not constructing this right...?

    uncodedLen = len(uncoded)

    kmpPartials = computeKMPPartials(uncoded)

    bestIndex = 0
    bestLength = 0
    m = 0
    i = 0

    while m < SEARCH_SIZE:
        if uncoded[i] == buffer[m + i + windowStart]:
            i +=1

            if(i == uncodedLen):
                bestLength = uncodedLen
                bestIndex = m + windowStart
                return bestIndex, bestLength
            
        else:
            if i > bestLength:
                bestLength = i
                bestIndex = m + windowStart
            
            m = m + i - kmpPartials[i]

            if kmpPartials[i] > 0:
                i = kmpPartials[i]
            else:
                i = 0

    return bestIndex, bestLength



def computeKMPPartials(pat):
    patlen = len(pat)
    kmpPartials = [0] * patlen
    kmpPartials[0] = -1

    i = 2
    j = 0

    while i < patlen:
        if pat[i-1] == pat[j]:
            # Candidate found. Continue matching to find longest partial
            j += 1
            kmpPartials[i] = j
            i += 1
        elif j > 0:
            # No more matches, but at least one candidate for fallback
            j = kmpPartials[j]
        else:
            # Never found a match
            kmpPartials[i] = 0
            i += 1

    return kmpPartials



# Open a file for reading and writing. If the file doesn't exist, create it.
# (Vanilla open() with mode "r+" raises an error if the file doesn't exist.)
def touchopen(filename, *args, **kwargs):
    fd = os.open(filename, os.O_RDWR | os.O_CREAT)
    return os.fdopen(fd, *args, **kwargs)



if __name__ == "__main__":

    # Check for incorrect usage.
    argc = len(sys.argv)
    if argc < 2 or argc > 4:
        print("Usage: {0:s} <inFile> [outFile] [outOffset]".format(
            sys.argv[0]
        ))
        sys.exit(1)

    # Copy the arguments.
    inFile = sys.argv[1]
    outFile = None
    if argc == 3 or argc == 4:
        outFile = sys.argv[2]
    outOffset = 0
    if argc == 4:
        outOffset = int(sys.argv[3], 16)

    # Read the input file.
    with open(inFile, "rb") as inStream:
        inBytes = bytearray(inStream.read())

    # Compress the data.
    outBytes = compress(inBytes)

    # Write the compressed output, if appropriate.
    if outFile is not None:
        with touchopen(outFile, "r+b") as outStream:
            outStream.seek(outOffset)
            outStream.write(outBytes)
            lastOffset = outStream.tell()
            print("Last offset written, inclusive: {0:X}".format(
                lastOffset - 1
            ))

    # Report statistics on the data.
    print("Uncompressed size: 0x{0:X} ({0:d}) bytes".format(len(inBytes)))
    print("Compressed size: 0x{0:X} ({0:d}) bytes".format(len(outBytes)))
    print("Ratio: {0:f}".format(len(outBytes) / len(inBytes)))

    # Exit.
    sys.exit(0)
