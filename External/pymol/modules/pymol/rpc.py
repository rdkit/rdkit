""" an XML-RPC server to allow remote control of PyMol
 
  Author: Greg Landrum (glandrum@users.sourceforge.net)
  Created:       January 2002
  $LastChangedDate$
  License:  This file is part of the RDKit. The contents are covered by the terms of
            the BSD license which is included in the file license.txt, found at the
            root of the RDKit source tree. 
  Requires:
            - a python xmlrpclib distribution containing the SimpleXMLRPCServer
              module (1.0 or greater should be fine)
            - python with threading enabled  
 
  RD Version: $Rev$            
"""

import SimpleXMLRPCServer
import threading, sys, time, types, os, tempfile
from pymol import cmd, cgo

# initial port to try for the server
_xmlPort = 9123
# number of alternate ports to try if the first fails
_nPortsToTry = 5


def rpcCmd(cmdText):
  """ executes a PyMol API command
 
   return value is either the result of the command or the empty string
 
  """
  res = cmd.do(cmdText)
  if res is not None:
    return res
  else:
    return ''


def rpcQuit():
  """ causes PyMol to quit """
  cmd.quit()
  return 1


def rpcZoom(what=''):
  """ executes cmd.zoom(what) """
  cmd.zoom(what)
  return 1


def rpcSet(prop, val, obj):
  """ executes a PyMol set command
 
   return value is either the result of the command or the empty string
 
  """
  res = cmd.set(prop, val, obj)
  if res is not None:
    return res
  else:
    return ''


def rpcGet(prop, obj):
  """ executes a PyMol get command
 
   return value is either the result of the command or the empty string
 
  """
  res = cmd.get(prop, obj)
  if res is not None:
    return res
  else:
    return ''


def rpcPing():
  """ Used to establish whether or not the server is alive.
 
  This is a good thing to call after establishing a connection just to
  make sure that everything is ok.
 
    Returns 1
 
  """
  return 1


def rpcLabel(pos, labelText, id='lab1', color=(1, 1, 1)):
  """ create a text label
 
    Arguments:
      pos: a 3 tuple with the position of the label
      text: a string with the label
      color: a 3 tuple with the color of the label. (1,1,1) is white
      id: (OPTIONAL) the name of the object to be created
 
    NOTE:
      at the moment this is, how you say, a hack
 
  """
  x, y, z = pos
  text = """
Atom
 
  1  0  0  0  0  0  0  0  0  0999 V2000
% 10.4f% 10.4f%10.4f C   0  0  0  0  0  0  0  0  0  0  0  0
M  END""" % (x, y, z)
  cmd.read_molstr(text, id)
  cmd.label("%s" % (id), '"%s"' % labelText)
  cmd.hide("nonbonded", id)
  cmd.set_color("%s-color" % id, color)
  cmd.color("%s-color" % id, id)
  return 1


def rpcResetCGO(id):
  """ removes a CGO from the local dictionary
 
  """
  global cgoDict
  if id == "*":
    cgoDict = {}
    res = 1
  elif id in cgoDict:
    del (cgoDict[id])
    res = 1
  else:
    res = 0
  return res


def rpcSphere(pos, rad, color, id='cgo', extend=1, transparent=0, transparency=0.5):
  """ create a sphere
 
    Arguments:
      pos: a 3 tuple with the position of the sphere
      rad: a float with the radius
      color: a 3 tuple with the color of the sphere. (1,1,1) is white
      id: (OPTIONAL) the name of the object to be created
      extend: (OPTIONAL) if this is nonzero, the object will be cleared
        before adding the new sphere.  Otherwise the sphere is appended
        to the ojbect
      transparent: (OPTIONAL) sets the object to be transparent
      transparency: (OPTIONAL) the percent transparency of the object
  """
  r, g, b = color
  x, y, z = pos
  if extend:
    obj = cgoDict.get(id, [])
  else:
    obj = []
  if not transparent:
    o = []
  else:
    o = [cgo.ALPHA, 1 - transparency]
  o.extend([cgo.COLOR, r, g, b, cgo.SPHERE, x, y, z, rad])
  obj.extend(o)
  cgoDict[id] = obj
  cmd.load_cgo(obj, id, 1)
  return 1


def rpcRenderCGO(cgoV, id='cgo', extend=1):
  """ renders a CGO vector
 
    Arguments:
      cgoV: a vector of floats
      id: (OPTIONAL) the name of the object to be created
      extend: (OPTIONAL) if this is nonzero, the object will be cleared
        before adding the new sphere.  Otherwise the sphere is appended
        to the ojbect
  """
  if extend:
    obj = cgoDict.get(id, [])
  else:
    obj = []
  obj.extend(cgoV)
  cmd.load_cgo(obj, id, 1)
  return 1


def rpcSpheres(sphereD, id='cgo', extend=1):
  """ create a sphere
 
    Arguments:
      sphereD: a series of (pos,rad,color,transparent,transparency) tuples
      id: (OPTIONAL) the name of the object to be created
      extend: (OPTIONAL) if this is nonzero, the object will be cleared
        before adding the new sphere.  Otherwise the sphere is appended
        to the ojbect
  """
  if extend:
    obj = cgoDict.get(id, [])
  else:
    obj = []
  for pos, rad, color, transparent, transparency in sphereD:
    r, g, b = color
    x, y, z = pos
    if not transparent:
      o = []
    else:
      o = [cgo.ALPHA, 1 - transparency]
    o.extend([cgo.COLOR, r, g, b, cgo.SPHERE, x, y, z, rad])
    obj.extend(o)
  cgoDict[id] = obj
  cmd.load_cgo(obj, id, 1)
  return 1


def rpcCylinder(end1, end2, rad, color1, id='cgo', color2=None, extend=1, transparent=0,
                transparency=0.5):
  """ create a cylinder
 
    Arguments:
      end1: a 3 tuple with the position of end1 of the sphere
      end2: a 3 tuple with the position of end1 of the sphere
      rad: a float with the radius
      color1: a 3 tuple with the color of end1 of the sphere. (1,1,1) is white
      id: (OPTIONAL) the name of the object to be created
      color2: (OPTIONAL) a 3 tuple with the color of end2 of the sphere. (1,1,1) 
is white
      extend: (OPTIONAL) if this is nonzero, the object will be cleared
        before adding the new sphere.  Otherwise the sphere is appended
        to the ojbect
      transparent: (OPTIONAL) sets the object to be transparent
      transparency: (OPTIONAL) the percent transparency of the object
 
    NOTE: the reason that color2 follows id is that I think clients are
    going to be interested in setting the id more often than they are going
    to care about the second color.
 
  """
  global cgoDict

  if color2 is None:
    color2 = color1
  r1, g1, b1 = color1
  r2, g2, b2 = color2
  x1, y1, z1 = end1
  x2, y2, z2 = end2
  if extend:
    obj = cgoDict.get(id, [])
  else:
    obj = []
  if not transparent:
    o = []
  else:
    o = [cgo.ALPHA, 1 - transparency]
  o.extend([cgo.CYLINDER,
            x1,
            y1,
            z1,
            x2,
            y2,
            z2,
            rad,
            r1,
            g1,
            b1,
            r2,
            g2,
            b2, ])
  obj.extend(o)
  cgoDict[id] = obj
  cmd.load_cgo(obj, id, 1)
  return 1


def rpcShow(objs):
  """ shows (enables) an object (or objects)"""
  if type(objs) not in (types.ListType, types.TupleType):
    objs = (objs, )

  for objName in objs:
    try:
      cmd.enable(objName)
    except Exception:
      res = 0
      break
    else:
      res = 1
  return res


def rpcHide(objs):
  """ hides (disables) an object (or objects) """
  if type(objs) not in (types.ListType, types.TupleType):
    objs = (objs, )

  for objName in objs:
    try:
      cmd.disable(objName)
    except Exception:
      res = 0
      break
    else:
      res = 1
  return res


def rpcDeleteObject(objName):
  """ deletes an object """
  try:
    cmd.delete(objName)
  except Exception:
    res = 0
  else:
    res = 1
  return res


def rpcDeleteAll():
  """ deletes all objects """
  res = cmd.delete('all')
  if res is not None:
    return res
  else:
    return ''


def colorObj(objName, colorScheme):
  """ sets an molecule's color scheme
    Arguments:
      - objName: the object (molecule) to change
      - colorScheme: name of the color scheme to use
        for the object (should be either 'std' or one of the
        color schemes defined in pymol.utils)
  """
  if colorScheme:
    if colorScheme == 'std':
      # this is an adaptation of the cbag scheme from util.py, but
      # with a gray carbon.
      cmd.color("magenta", "(" + objName + ")", quiet=1)
      cmd.color("oxygen", "(elem O and " + objName + ")", quiet=1)
      cmd.color("nitrogen", "(elem N and " + objName + ")", quiet=1)
      cmd.color("sulfur", "(elem S and " + objName + ")", quiet=1)
      cmd.color("hydrogen", "(elem H and " + objName + ")", quiet=1)
      cmd.color("gray", "(elem C and " + objName + ")", quiet=1)
    elif hasattr(utils, colorScheme):
      fn = getattr(utils, colorScheme)
      fn(objName, quiet=1)
    res = 1
  else:
    res = 0
  return res


def rpcLoadPDB(data, objName, colorScheme='', replace=1):
  """ loads a molecule from a pdb string
 
    Arguments:
      data: the mol block
      objName: name of the object to create
      colorScheme: (OPTIONAL) name of the color scheme to use
        for the molecule (should be either 'std' or one of the
        color schemes defined in pymol.utils)
      replace: (OPTIONAL) if an object with the same name already
        exists, delete it before adding this one
 
 
  """
  from pymol import util
  if replace:
    cmd.delete(objName)
  res = cmd.read_pdbstr(data, objName)
  colorObj(objName, colorScheme)
  if res is not None:
    return res
  else:
    return ''


def rpcLoadMolBlock(data, objName, colorScheme='', replace=1):
  """ loads a molecule from a mol block
 
    Arguments:
      data: the mol block
      objName: name of the object to create
      colorScheme: (OPTIONAL) name of the color scheme to use
        for the molecule (should be either 'std' or one of the
        color schemes defined in pymol.utils)
      replace: (OPTIONAL) if an object with the same name already
        exists, delete it before adding this one
 
  """
  from pymol import util
  if replace:
    cmd.delete(objName)
  res = cmd.read_molstr(data, objName)
  colorObj(objName, colorScheme)
  if res is not None:
    return res
  else:
    return ''


def rpcLoadFile(fileName, objName='', format='', colorScheme='', replace=1):
  """ loads an object from a file
 
    Arguments:
      fileName: the file to load
      objName: (OPTIONAL) name of the object to create
      format: (OPTIONAL) the format of the input file
      colorScheme: (OPTIONAL) name of the color scheme to use
        for the object (should be either 'std' or one of the
        color schemes defined in pymol.utils)
      replace: (OPTIONAL) if an object with the same name already
        exists, delete it before adding this one
  """
  if not objName:
    objName = fileName.split('.')[0]
  if replace:
    cmd.delete(objName)
  res = cmd.load(fileName, objName, format=format)
  colorObj(objName, colorScheme)
  if res is not None:
    return res
  else:
    return ''


def rpcLoadSurface(fileName, objName, format='', surfaceLevel=1.0):
  """ loads surface data from a file and adds an isosurface
 
    Arguments:
      fileName: the file to load
      objName: (OPTIONAL) name of the object to create
      format: (OPTIONAL) the format of the input file
      surfaceLevel: (OPTIONAL) the isosurface level
 
   """
  if not objName:
    objName = fileName.split('.')[0]
  gridName = 'grid-%s' % objName
  res = cmd.load(fileName, gridName, format='')
  cmd.isosurface(objName, gridName, level=surfaceLevel)
  if res is not None:
    return res
  else:
    return ''


def rpcLoadSurfaceData(data, objName='surface', format='', surfaceLevel=1.0):
  """ loads surface data from a string and adds an isosurface
 
    Arguments:
      data: the data to load
      objName: (OPTIONAL) name of the object to create
      format: (OPTIONAL) the format of the input file
      surfaceLevel: (OPTIONAL) the isosurface level
 
   """
  gridName = 'grid-%s' % objName
  # it would be nice if we didn't have to go by way of the temporary file,
  # but at the moment pymol will only read shapes from files
  with tempfile.NamedTemporaryFile('w+', suffix='.grd', delete=False) as tmp:
    tmp.write(data)
  res = rpcLoadSurface(tmp.name, objName, format='', surfaceLevel=surfaceLevel)
  os.unlink(tmp.name)
  if res is not None:
    return res
  else:
    return ''


def rpcSave(filename, objName='all', state=0, format=''):
  """ executes a cmd.save command
 
    Arguments:
     - filename: output filename
     - objName: (OPTIONAL) object(s) to be saved
     - state: (OPTIONAL)
     - format: (OPTIONAL) output format
 
  """
  res = cmd.save(filename, objName, state, format)
  if res is not None:
    return res
  else:
    return ''


def rpcRotate(vect, objName='', state=-1):
  """ rotates objects
 
    Arguments:
     - vect: a sequence with x y and z rotations
     - objName: (OPTIONAL) object to be rotated
     - state: (OPTIONAL) if zero only visible states are rotated,
       if -1 (the default), all states are rotated
 
  """
  cmd.rotate('x', vect[0], objName, state=state)
  cmd.rotate('y', vect[1], objName, state=state)
  cmd.rotate('z', vect[2], objName, state=state)
  return 1


def rpcTranslate(vect, objName='all', state=-1):
  """ translates objects
 
    Arguments:
     - vect: a sequence with x y and z translations
     - objName: (OPTIONAL) object to be translated
     - state: (OPTIONAL) if zero only visible states are translated,
       if -1 (the default), all states are translated
  """
  cmd.translate(vect, objNAme, state=state)
  return 1


def rpcGetNames(what='selections', enabledOnly=1):
  """ returns the results of cmd.get_names(what) """
  return cmd.get_names(what, enabled_only=enabledOnly)


def rpcIdentify(what='all', mode=0):
  """ returns the results of cmd.identify(what,mode) """
  return cmd.identify(what, mode=mode)


def rpcIndex(what='all'):
  """ returns the results of cmd.index(what) """
  return cmd.index(what)


def rpcCountAtoms(what='all'):
  """ returns the results of cmd.count_atoms(what) """
  return cmd.count_atoms(what)


def rpcIdAtom(what='all', mode=0):
  """ returns the results of cmd.id_atom(what) """
  return cmd.id_atom(what, mode=mode)


def rpcGetAtomCoords(what='all', state=0):
  """ returns the results of cmd.get_atom_coords(what,state) """
  return cmd.get_atom_coords(what, state=state)


def rpcHelp(what=''):
  """ returns general help text or help on a particular command """
  global serv
  res = 'Command Not Found'
  if not what:
    res = serv.funcs.keys()
  else:
    funcs = serv.funcs
    if what in funcs:
      fn = funcs[what]
      res = "Function: %s(" % what
      defs = fn.func_defaults
      if defs:
        code = fn.func_code
        nDefs = len(defs)
        args = []
        i = -1
        for i in range(code.co_argcount - nDefs):
          args.append(code.co_varnames[i])
        for j in range(nDefs):
          vName = code.co_varnames[j + i + 1]
          args.append("%s=%s" % (vName, repr(defs[j])))
        res += ','.join(args)
      res += ')\n'
      if fn.func_doc:
        res += fn.func_doc
  return res


def launch_XMLRPC(hostname='', port=_xmlPort, nToTry=_nPortsToTry):
  """ launches the xmlrpc server into a separate thread
 
    Arguments:
      hostname: (OPTIONAL) name of the host for the server
         (defaults to be the name of the localhost)
      port: (OPTIONAL) the first port to try for the server
      nToTry: (OPTIONAL) the number of possible ports to try
        (in case the first can't be opened)
 
  """
  if not hostname:
    import os
    hostname = os.environ.get('PYMOL_RPCHOST', '')
    if not hostname or hostname.upper() == 'LOCALHOST':
      hostname = 'localhost'
    else:
      import socket
      hostname = socket.gethostbyname(socket.gethostname())

  global cgoDict, serv
  cgoDict = {}
  for i in range(nToTry):
    try:
      serv = SimpleXMLRPCServer.SimpleXMLRPCServer((hostname, port + i), logRequests=0)
    except Exception:
      serv = None
    else:
      break
  if serv:
    print('xml-rpc server running on host %s, port %d' % (hostname, port + i))
    serv.register_function(rpcCmd, 'do')
    serv.register_function(rpcQuit, 'quit')
    serv.register_function(rpcSet, 'set')
    serv.register_function(rpcGet, 'get')
    serv.register_function(rpcPing, 'ping')
    serv.register_function(rpcResetCGO, 'resetCGO')
    serv.register_function(rpcRenderCGO, 'renderCGO')
    serv.register_function(rpcSphere, 'sphere')
    serv.register_function(rpcSpheres, 'spheres')
    serv.register_function(rpcCylinder, 'cylinder')
    serv.register_function(rpcHide, 'hide')
    serv.register_function(rpcShow, 'show')
    serv.register_function(rpcZoom, 'zoom')
    serv.register_function(rpcDeleteObject, 'deleteObject')
    serv.register_function(rpcDeleteAll, 'deleteAll')
    serv.register_function(rpcLoadPDB, 'loadPDB')
    serv.register_function(rpcLoadMolBlock, 'loadMolBlock')
    serv.register_function(rpcLoadSurface, 'loadSurface')
    serv.register_function(rpcLoadSurfaceData, 'loadSurfaceData')
    serv.register_function(rpcLoadFile, 'loadFile')
    serv.register_function(rpcSave, 'save')
    serv.register_function(rpcLabel, 'label')
    serv.register_function(rpcRotate, 'rotate')
    serv.register_function(rpcTranslate, 'translate')
    serv.register_function(rpcGetNames, 'getNames')
    serv.register_function(rpcIdentify, 'identify')
    serv.register_function(rpcIndex, 'index')
    serv.register_function(rpcCountAtoms, 'countAtoms')
    serv.register_function(rpcIdAtom, 'idAtom')
    serv.register_function(rpcHelp, 'help')
    serv.register_function(rpcGetAtomCoords, 'getAtomCoords')
    serv.register_introspection_functions()
    t = threading.Thread(target=serv.serve_forever)
    t.setDaemon(1)
    t.start()
  else:
    print('xml-rpc server could not be started')
