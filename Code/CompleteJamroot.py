import os,sys
try:
    import numpy
except ImportError:
    print >>sys.stderr,"ERROR: numpy module not found.\nThe RDKit python wrappers will not build correctly, but the rest of the code should be fine."
    sys.exit(0)

incDir = os.path.join(os.path.split(numpy.__file__)[0],'core','include')
inD = file('Jamroot.in','r').read()
print >>file('Jamroot','w+'),inD.replace('NUMPY_INCLUDE_DIR',incDir)
    
