"""
Script provides a specialized edit distance for partial cognates.

Script requires: lingpy (http://lingpy.org) to calculate the edit distance
between strings.
"""
from lingpy import *

def pcdist(seqA, seqB, insert=6, delete=3):
    """
    Calculate a specified distance between two strings.
    """

    almA, almB, score = nw_align(seqA, seqB)
    
    score = 0
    matches = 0

    if seqA == '-' and seqB == '-':
        return 0
    elif seqA == '-':
        return insert
    elif seqB == '-':
        return delete
    
    for a,b in zip(almA, almB):
        if a == '-':
            score += 1
        elif b == '-':
            score += 2
        elif a != b:
            score += 2
        elif a == b:
            matches += 1

    return score

def editdist(seqA, seqB, insert=6, delete=3):
    
    score = 0


    d = edit_dist(seqA, seqB)

    if seqA == '-' and seqB != '-':
        return insert
    elif seqA != '-' and seqB == '-':
        return delete
    
    almA, almB, x = nw_align(seqA, seqB)

    for a,b in zip(almA, almB):
        if a == b:
            pass
        elif a == '-' or b == '-':
            score += 1
        else:
            score += 2
    
    return score

def simple(seqA, seqB, insert=6, delete=3):

    if seqA == '-' and seqB != '-':
        return insert
    elif seqA != '-' and seqB == '-':
        return delete

    if seqA != seqB:
        return 1
    else:
        return 0


if __name__ == '__main__':

    patterns = [
        "背脊",
        "背脢",
        "背心",
        "背",
        "背梁",
        "脊背",
        "背脊",
        "月聘",
        "脊梁",
        "背脊",
        "背脊",
        "背瑯",
        "背",
        "目脊",
        "胶脊",
        "背脊",
        "背脊",
        "胛脊",
        "背",
        "脊梁",
        "背脊身",
        "背",
        "背心",
        "脊背",
        "脊背",
        "背",
        "背心",
        "背",
        "背心"
        ]
    patterns = sorted(set(patterns))
    for i,p1 in enumerate(patterns):
        for j,p2 in enumerate(patterns):
            if i < j:
                s,m = pcdist(p1, p2)
                print(p1, '\t->\t', p2, '\t', s, '\t', m)


