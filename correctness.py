#!/usr/bin/env python

with open('statistics.csv', 'r') as f:
    buffered = f.read()

buffered = buffered.split('\n')
signatures = list()
for signature in buffered:
    signatures.append(signature.split('\t'))

for signature in signatures:
    if signature[0] == 'Smiles':
        continue
    if (float(signature[3]) + float(signature[4])) != (float(signature[5]) + float(signature[6])) and\
        (float(signature[3]) + float(signature[4])) != 0 and \
        (float(signature[5]) + float(signature[6]) != 0):
        print signature[0]
