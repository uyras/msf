#!/usr/bin/env python3

import numpy as np
import sys

if len(sys.argv) != 3:
    print(f"Использование: {sys.argv[0]} <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

def sq(arr):
    return arr.reshape((-1,int(np.sqrt(arr.shape[0]))))

pre = "this file is translated from long to short format\n"
with open(input_file) as f:
    for line in f:
        if line[0]=="#":
            pre += line[1:]
        else:
            break
pre = pre.strip()
dat = np.loadtxt(input_file)
rea = sq(dat[:,3])
np.savetxt(output_file,rea,header=pre,fmt='%.8g')