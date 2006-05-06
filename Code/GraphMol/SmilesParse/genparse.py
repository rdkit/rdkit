import os,sys

bindir = 'c:\\cygwin\\bin\\';
os.system('%sbison -p yy%s_ -t -d -o %s.tab.cpp %s.y' % (bindir,sys.argv[1],sys.argv[1],sys.argv[1]) )

