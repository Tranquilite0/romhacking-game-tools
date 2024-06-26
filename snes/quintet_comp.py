#!/usr/bin/env python3
#
# Quintet Compressor
# Osteoclave
# 2012-02-04
#
# The compression format is described in the decompressor.
#
# Modified by Tranquilite
# Corrected issue where algorithm would only match up to 15 bytes instead of
# the maximum possible of 17 (should match or beat the existing compression).
# Switched to KMP search for a little speed boost.
#
# This code uses python-bitstring:
# https://pypi.org/project/bitstring/

import os
import sys
import bitstring


# Define some useful constants.
SEARCH_LOG2 = 8
SEARCH_SIZE = 2 ** SEARCH_LOG2
LOOKAHEAD_LOG2 = 4
LOOKAHEAD_SIZE = 2 ** LOOKAHEAD_LOG2
BIT_PASTCOPY = 0
BIT_LITERAL = 1
MAX_UNCODED = 1
MAX_CODED = LOOKAHEAD_SIZE + MAX_UNCODED
MIN_CODED = 2



def compress(inBytes):
    # Prepare the memory buffer.
    inBuffer = bytearray(SEARCH_SIZE + len(inBytes))
    inBuffer[:SEARCH_SIZE] = [0x20] * SEARCH_SIZE
    inBuffer[SEARCH_SIZE:] = inBytes

    # Prepare for compression.
    output = bitstring.BitArray()
    output += bitstring.pack("uintle:16", len(inBytes))
    currentIndex = SEARCH_SIZE

    # Main compression loop.
    while currentIndex < len(inBuffer):
        bestIndex = 0
        bestLength = 0

        #bestIndex, bestLength = BruteForceSearch(inBuffer, currentIndex)
        bestIndex, bestLength = KMPSearch(inBuffer, currentIndex-SEARCH_SIZE, currentIndex )

        # Write the next block of compressed output.
        if bestLength >= 2:
            # For some reason, the decompressor expects the pastcopy
            # source values to be offset by 0xEF. I have no idea why.
            bestIndex = (bestIndex + 0xEF) & 0xFF
            output += bitstring.pack("uint:1", BIT_PASTCOPY)
            output += bitstring.pack(
                "uint:n=v", n = SEARCH_LOG2, v = bestIndex
            )
            output += bitstring.pack(
                "uint:n=v", n = LOOKAHEAD_LOG2, v = bestLength - 2
            )
            currentIndex += bestLength
        else:
            output += bitstring.pack("uint:1", BIT_LITERAL)
            output += bitstring.pack("uint:8", inBuffer[currentIndex])
            currentIndex += 1

    # Return the compressed data.
    return output.tobytes()


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
