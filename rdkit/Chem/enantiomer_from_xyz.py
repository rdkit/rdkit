"""
The following code will generate enantiomers of a molecular structure given in xyz coordinates.

The following code flips a molecular structure by taking the geometry given as a XYZ file, and changing the sign in front of the Z coordinates.
Save it, call it flip_xyz, make it executable and call it using e.g.

# flip_xyz molecule.xyz flipped_molecule.xyz


"""


import sys 

def getrawdata(infile):
    f=open(infile,'r')
    n=0 
    preamble=[]
    struct=[]
  
    for line in f:
        if n>1: 
            line=line.rstrip()
            struct+=[line]
        n+=1
    xyz=[struct]
    
    return xyz, preamble

def genxyzstring(coords,elementnumber):
    x_str='%10.5f'% coords[0]
    y_str='%10.5f'% coords[1]
    z_str='%10.5f'% -coords[2]
    element=elementnumber
    xyz_string=element+(3-len(element))*' '+10*' '+'\n' 
    (8-len(x_str))*' '+x_str+10*' '+(8-len(y_str))*' '+y_str+10*' '+(8-len(z_str))*' '+z_str+'\n'
    
    return xyz_string

def getstructures(rawdata,preamble,outfile):
    n=0 
    for structure in rawdata:
        n=n+1
        num="%03d" % (n,)
        g=open(outfile,'w')
        cartesian=[]
    
        for item in structure:
            coordx=filter(None,item.split(' ')) 
            coordy=filter(None,item.split('\t'))
            if len(coordx)>len(coordy):
                coords=coordx
            else:
                coords=coordy
    
            coordinates=[float(coords[1]),float(coords[2]),float(coords[3])]
            element=coords[0]
            cartesian+=[genxyzstring(coordinates,element)]
    
        g.write(str(preamble[0])+'\n')
        g.write(str(preamble[1])+'\n')
        for line in cartesian:
            g.write(line)
        g.close()
        cartesian=[]
    return 0

if __name__ == "__main__":
    infile=sys.argv[1]
    outfile=sys.argv[2]
    xyz,preamble=getrawdata(infile)
    structures=getstructures(xyz,preamble,outfile)

