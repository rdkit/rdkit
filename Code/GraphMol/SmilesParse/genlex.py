import os,sys

bindir = 'c:/cygwin/bin/';

if len(sys.argv) > 2:
    runDir = sys.argv[2]
    os.chdir(runDir)
os.system('%sflex -Pyy%s_ -olex.yy.tmp.cpp %s.ll' % (bindir,sys.argv[1],sys.argv[1] ))
os.system('%sgrep -v unistd lex.yy.tmp.cpp > lex.yy%s.cpp' % (bindir,sys.argv[1] ))
