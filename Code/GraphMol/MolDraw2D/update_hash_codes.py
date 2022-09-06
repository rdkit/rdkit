#!/usr/bin/env python

# Script to update the hash codes in the cpp codes.
# Takes 2 arguments, the output from running one of the test programs
# giving the values to update, and the corresponding cpp
# file.  Output is to {file.cpp}.update.
# Note that the PNG hash codes are only for Freetype builds, so it's
# better to do the non-Freetype builds first, otherwise the PNG
# hashes will definitely be wrong.

import json
import re
import sys


updates = {}
with open(sys.argv[1], 'r') as f:
    for line in f.readlines():
        if line.startswith('file'):
            line_bits = line.strip().split()
            updates[(line_bits[1], line_bits[8])] = (line_bits[1], line_bits[4])

re_str = r'{"([\w.-]*.\w*)", (\d*U)},'
patt = re.compile(re_str)

with open(sys.argv[2], 'r') as f:
    cpp_file = f.read()

hashes = patt.findall(cpp_file)

def hash_replace(match):
    match_bits = match.group().split(',')
    match_tuple = (match_bits[0][2:-1], match_bits[1][1:-1])
    if match_tuple in updates:
        rep_tuple = updates[match_tuple]
        rep_str = f'{{"{rep_tuple[0]}", {rep_tuple[1]}}},'
        print(f'replacing {match.group()} with {rep_str}')
        return rep_str
    else:
        return match.group()
    
for hash in hashes:
    print(hash)
    if hash in updates:
        print(hash, updates[hash])
    
new_cpp_file = patt.sub(hash_replace, cpp_file)

with open(f'{sys.argv[2]}.update', 'w') as f:
    f.write(new_cpp_file)
    
