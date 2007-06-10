import os,sys

if len(sys.argv) > 2:
    runDir = sys.argv[2]
    os.chdir(runDir)
    
bindir = 'c:\\cygwin\\bin\\';
os.system('%sbison -p yy%s_ -t -d -o %s.tab.cpp %s.y' % (bindir,sys.argv[1],sys.argv[1],sys.argv[1]) )

