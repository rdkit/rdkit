import sys, os

dicfile = open(sys.argv[1], 'r')
line = dicfile.readline()

start = 0
dict = {}
while line :
    tmp = line.strip()
    if (tmp == "#if 1") : 
        start = 1

    if (tmp == "#else") :
        break
    
    if (start) :
        tlst = tmp.split()
        if len(tlst) == 3 :
            dict[tlst[1]] = tlst[2]
    line = dicfile.readline()

dicfile.close()

filelst = open(sys.argv[2], 'r')

fname = filelst.readline().strip()

while fname :
    tmpfi = open(fname, 'r')
    oname = fname + ".tmp"
    tmpfo = open(oname, 'w')
    print "Working on : ", fname
    line = tmpfi.readline()
    while line :
        tmp = line.strip()
        tlst = tmp.split()
        if (len(tlst) == 2) :
            if (tlst[0] == "#include") :
                if (dict.has_key(tlst[1])) :
                    oline = tlst[0] + " " + dict[tlst[1]] + "\n"
                    line  = "//" + tmp + " // changed of VC++\n" + oline
        tmpfo.write(line)
        line = tmpfi.readline()
    
    tmpfi.close()
    tmpfo.close()

    os.system("mv " + oname + " " + fname)
    fname = filelst.readline().strip()

filelst.close()
