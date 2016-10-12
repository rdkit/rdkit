#!/usr/bin/env python
import sys

for file in sys.argv[1:]:
    if os.path.exists(file):
        text = open(file).read()
        text = text.replace("|| (LONG_MAX == INT_MAX)", "")
        open(file,'w').write(text)
