
bison -t -d smiles.y 
mv smiles.tab.c smiles.tab.cpp
flex smiles.flex
grep -v unistd lex.yy.c | sed s/lex.yy.c/lex.yy.cpp/ > lex.yy.cpp
