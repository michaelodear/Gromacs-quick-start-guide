#!/usr/bin/python

import numpy

# Ask for residue of interest

res = int(input("What residue do you want to look at?: "))

test = input("Name of test file: ")
control = input("Name of control to compare against: ")
output = input("Name of output plot: ")

# Create 2D arrays and initialise to 0
residues_WT = numpy.zeros([11,8])
residues_compound = numpy.zeros([11,8])

def parse_data(residues, file):
    with open(file) as f:
        i = 0
        for line in f:
            check = (i - res) % 121
            if (check < 11):
                line = line.split()
                for j in range(8):
                    residues[check][j] += 100*float(line[j+1])
            i += 1

    return residues

# Fill arrays
WT_structures = parse_data(residues_WT, control)
compound_structures = parse_data(residues_compound, test)

# Average and difference

difference = []

for x, y in zip(WT_structures, compound_structures):
    diff = []
    for m,n in zip(x,y):
        diff.append((n-m)/14)
    difference.append(diff)

a = open(output, 'w+')

i = 0

header = open(control, "r")

for line in header:
    if(i<13):
        a.write(line)
    else:
        break
    i += 1

k = -5

for entry in difference:
    test = ''.join(str(entry))
    test = test[1:-1]
    test = test.replace(',','')
    a.write(str(k)+' '+test)
    a.write('\n')
    k += 1