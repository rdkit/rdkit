#
#  Copyright (C) 2002  Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#

""" undo functionality

"""

class UndoError(Exception):
  def __init__(self,value):
    self.value = value
  def __str__(self):
    return `self.value`


class UndoObj:
  def __init__(self,name,undoFunc,undoArgs,redoFunc,redoArgs):
    self._name = name
    self._uFunc = undoFunc
    self._uArgs = undoArgs
    self._rFunc = redoFunc
    self._rArgs = redoArgs
  def undo(self):
    #print 'UNDO:',self.name()
    apply(self._uFunc,self._uArgs)
  def redo(self):
    #print 'REDO:',self.name()
    apply(self._rFunc,self._rArgs)
  def accepts(self,other):
    return 0
  def accept(self,other):
    pass
  def name(self):
    return self._name

class UndoStack:
  """

    self._pos points at the thing to be undone, for redo look forward

  """
  def __init__(self,size=100):
    self._size = size
    self._stack = [None]*self._size
    self._pos = -1

  def dbg(self):
    for entry in self._stack:
      if entry: print '"%s"'%(entry._name),
      else: print '*',
    print  
  def capacity(self):
    return self._size
  
  def undo(self):
    if self._pos < 0 or self._stack[self._pos] is None:
      raise UndoError,'empty stack'
    else:
      self._stack[self._pos].undo()
      self._pos -= 1
      return 1

  def undoQuery(self):
    if self._pos < 0 or self._stack[self._pos] is None:
      res = None
    else:
      res = self._stack[self._pos].name()
    return res  

  def redo(self):
    self._pos += 1
    if self._pos > self._size-1 or self._stack[self._pos] is None:
      self._pos -= 1
      raise UndoError,'empty stack'
    else:
      self._stack[self._pos].redo()
      return 1

  def redoQuery(self):
    self._pos += 1
    if self._pos > self._size-1 or self._stack[self._pos] is None:
      res = None
    else:
      res = self._stack[self._pos].name()
    self._pos -= 1
    return res

  def push(self,obj):
    # see if we can just glom on the request
    if self._stack[self._pos] and self._stack[self._pos].accepts(obj):
      self._stack[self._pos].accept(obj)
    else:
      self._pos += 1
      if self._pos == self._size:
        del self._stack[0]
        self._stack.append(None)
        self._pos -= 1
      self._stack[self._pos] = obj
      for i in range(self._pos+1,self._size):
        self._stack[i] = None

  def reset(self):
    self._pos = 0

if __name__ == '__main__':
  import unittest
  class TestCase(unittest.TestCase):
    def setUp(self):
      self.s = UndoStack(4)
      self.var = [1]
    def _fn1(self,x=None):
      if x is None: x = self.var
      x[0] += 1
    def _fn2(self,x=None):
      if x is None: x = self.var
      x[0] -= 1
    def test1(self):
      self.s = UndoStack()
      try:
        self.s.undo()
      except UndoError,msg:
        ok = 1
      else:
        ok = 0
      assert ok,'failed to generate expected exception'  
    def test2(self):
      self.s = UndoStack(size=4)
      self.s.push(UndoObj('test',self._fn2,(self.var,),self._fn1,(self.var,)))
      self._fn1()
      assert self.var==[2],'bad starting conditions'
      self.s.undo()
      assert self.var==[1],'bad undo result: %s'%(self.var)
      self.s.redo()
      assert self.var==[2],'bad redo result'
      try:
        self.s.redo()
      except UndoError,msg:
        ok = 1
      else:
        ok = 0
      assert ok,'redo failed to generate expected exception'  
    def test3(self):
      """ multiple undos, redos """
      self.s = UndoStack(size=4)
      expect = self.var[:]
      for i in range(3):
        self.s.push(UndoObj('test',self._fn2,(self.var,),self._fn1,(self.var,)))
        self._fn1()
        self._fn1(expect)
        assert expect==self.var,'bad start match'
      for i in range(3):
        self.s.undo()
        self._fn2(expect)
        assert expect==self.var,'bad start match'
      for i in range(3):
        self.s.redo()
        self._fn1(expect)
        assert expect==self.var,'bad start match'
    def test4(self):
      """  stack rolls """
      self.s = UndoStack(size=4)
      expect = self.var[:]
      for i in range(10):
        self.s.push(UndoObj('test',self._fn2,(self.var,),self._fn1,(self.var,)))
        self._fn1()
        self._fn1(expect)
        assert expect==self.var,'bad start match'
      for i in range(3):
        self.s.undo()
        self._fn2(expect)
        assert expect==self.var,'bad start match'
      for i in range(3):
        self.s.redo()
        self._fn1(expect)
        assert expect==self.var,'bad start match'
      
    def test5(self):
      """  queries and stack rolls """
      self.s = UndoStack(size=4)
      expect = self.var[:]
      for i in range(10):
        self.s.push(UndoObj('%d'%i,self._fn2,(self.var,),self._fn1,(self.var,)))
        self._fn1()
        assert self.s.undoQuery()=='%d'%i,'bad undo query: %s'%self.s.undoQuery()
      pos = i
      for i in range(3):
        assert self.s.undoQuery()=='%d'%pos,'bad undo query: %s'%self.s.undoQuery()
        self.s.undo()
        pos -= 1
      pos += 1
      for i in range(3):
        assert self.s.redoQuery()=='%d'%pos,'bad redo query: %s'%self.s.redoQuery()
        self.s.redo()
        pos += 1

      
    def test6(self):
      """  partial rollback followed by revision """
      self.s = UndoStack(size=4)
      expect = self.var[:]
      for i in range(10):
        self.s.push(UndoObj('%d'%i,self._fn2,(self.var,),self._fn1,(self.var,)))
        self._fn1()
        assert self.s.undoQuery()=='%d'%i,'bad undo query: %s'%self.s.undoQuery()
      pos = i
      for i in range(3):
        assert self.s.undoQuery()=='%d'%pos,'bad undo query: %s'%self.s.undoQuery()
        self.s.undo()
        pos -= 1

      for i in range(3):
        self.s.push(UndoObj('%d'%i,self._fn2,(self.var,),self._fn1,(self.var,)))
        self._fn1()
        assert self.s.undoQuery()=='%d'%i,'bad undo query: %s'%self.s.undoQuery()

    # FIX: test accepts/accept

      
  unittest.main()

