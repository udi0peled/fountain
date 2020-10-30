#!/usr/bin/python3
import math

def minimal_val_for_bitlen(bitlen):
    curr = 0
    for i in reversed(range(bitlen)):
        curr |= (1 << i)
        if (curr * math.floor(math.log(curr,2)) > bitlen):
            curr &= ~(1 << i)
    return (curr, math.floor(math.log(curr, 2)))
