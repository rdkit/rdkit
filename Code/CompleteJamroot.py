import os,sys
inD = file('Jamroot.in','r').read()
try:
    import numpy
except ImportError:
    print >>sys.stderr,"ERROR: numpy module not found.\nThe RDKit python wrappers will not build correctly, but the rest of the code should be fine."
else:
    numpyDir=os.path.split(numpy.__file__)[0].replace('\\','/')
    incDir = '/'.join((numpyDir,'core','include'))
    inD = inD.replace('NUMPY_INCLUDE_DIR',incDir)
print >>file('Jamroot','w+'),inD
    
