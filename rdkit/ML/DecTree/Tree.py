#
#  Copyright (C) 2000-2008  greg Landrum and Rational Discovery LLC
#
""" Implements a class used to represent N-ary trees

"""


import pickle


# FIX: the TreeNode class has not been updated to new-style classes
# (RD Issue380) because that would break all of our legacy pickled
# data. Until a solution is found for this breakage, an update is
# impossible.
class TreeNode:
  """ This is your bog standard Tree class.

   the root of the tree is just a TreeNode like all other members.
  """

  def __init__(self, parent, name, label=None, data=None, level=0, isTerminal=0):
    """ constructor

     **Arguments**

       - parent: the parent of this node in the tree

       - name: the name of the node

       - label: the node's label (should be an integer)

       - data: an optional data field

       - level: an integer indicating the level of this node in the hierarchy
         (used for printing)

       - isTerminal: flags a node as being terminal.  This is useful for those
         times when it's useful to know such things.

    """
    self.children = []
    self.parent = parent
    self.name = name
    self.data = data
    self.terminalNode = isTerminal
    self.label = label
    self.level = level
    self.examples = []

  def NameTree(self, varNames):
    """ Set the names of each node in the tree from a list of variable names.

     **Arguments**

       - varNames: a list of names to be assigned

     **Notes**

        1) this works its magic by recursively traversing all children

        2) The assumption is made here that the varNames list can be indexed
           by the labels of tree nodes

    """
    if self.GetTerminal():
      return
    else:
      for child in self.GetChildren():
        child.NameTree(varNames)
      self.SetName(varNames[self.GetLabel()])

  NameModel = NameTree

  def AddChildNode(self, node):
    """ Adds a TreeNode to the local list of children

     **Arguments**

       - node: the node to be added

     **Note**

       the level of the node (used in printing) is set as well

    """
    node.SetLevel(self.level + 1)
    self.children.append(node)

  def AddChild(self, name, label=None, data=None, isTerminal=0):
    """ Creates a new TreeNode and adds a child to the tree

      **Arguments**

       - name: the name of the new node

       - label: the label of the new node (should be an integer)

       - data: the data to be stored in the new node

       - isTerminal: a toggle to indicate whether or not the new node is
         a terminal (leaf) node.

      **Returns*

        the _TreeNode_ which is constructed

    """
    child = TreeNode(self, name, label, data, level=self.level + 1, isTerminal=isTerminal)
    self.children.append(child)
    return child

  def PruneChild(self, child):
    """ Removes the child node

      **Arguments**

        - child: a TreeNode

    """
    self.children.remove(child)

  def ReplaceChildIndex(self, index, newChild):
    """ Replaces a given child with a new one

      **Arguments**

        - index: an integer

        - child: a TreeNode

    """
    self.children[index] = newChild

  def GetChildren(self):
    """ Returns a python list of the children of this node

    """
    return self.children

  def Destroy(self):
    """ Destroys this node and all of its children

    """
    for child in self.children:
      child.Destroy()
    self.children = []
    # clean up circular references
    self.parent = None

  def GetName(self):
    """ Returns the name of this node

    """
    return self.name

  def SetName(self, name):
    """ Sets the name of this node

    """
    self.name = name

  def GetData(self):
    """ Returns the data stored at this node

    """
    return self.data

  def SetData(self, data):
    """ Sets the data stored at this node

    """
    self.data = data

  def GetTerminal(self):
    """ Returns whether or not this node is terminal

    """
    return self.terminalNode

  def SetTerminal(self, isTerminal):
    """ Sets whether or not this node is terminal

    """
    self.terminalNode = isTerminal

  def GetLabel(self):
    """ Returns the label of this node

    """
    return self.label

  def SetLabel(self, label):
    """ Sets the label of this node (should be an integer)

    """
    self.label = label

  def GetLevel(self):
    """ Returns the level of this node

    """
    return self.level

  def SetLevel(self, level):
    """ Sets the level of this node

    """
    self.level = level

  def GetParent(self):
    """ Returns the parent of this node

    """
    return self.parent

  def SetParent(self, parent):
    """ Sets the parent of this node

    """
    self.parent = parent

  def Print(self, level=0, showData=0):
    """ Pretty prints the tree

      **Arguments**

        - level: sets the number of spaces to be added at the beginning of the output

        - showData: if this is nonzero, the node's _data_ value will be printed as well

      **Note**

        this works recursively

    """
    if showData:
      print('%s%s: %s' % ('  ' * level, self.name, str(self.data)))
    else:
      print('%s%s' % ('  ' * level, self.name))

    for child in self.children:
      child.Print(level + 1, showData=showData)

  def Pickle(self, fileName='foo.pkl'):
    """ Pickles the tree and writes it to disk

    """
    with open(fileName, 'wb+') as pFile:
      pickle.dump(self, pFile)

  def __str__(self):
    """ returns a string representation of the tree

      **Note**

        this works recursively

    """
    here = '%s%s\n' % ('  ' * self.level, self.name)
    for child in self.children:
      here = here + str(child)
    return here

  def __cmp__(self, other):
    """ allows tree1 == tree2

      **Note**

        This works recursively
    """
    return (self < other) * -1 or (other < self) * 1

  def __lt__(self, other):
    """ allows tree1 < tree2

      **Note**

        This works recursively
    """
    try:
      nChildren = len(self.children)
      oChildren = len(other.children)
      if str(type(self)) < str(type(other)):
        return True
      if self.name < other.name:
        return True
      if self.label is not None:
        if other.label is not None:
          if self.label < other.label:
            return True
        else:
          return False
      elif other.label is not None:
        return True
      if nChildren < oChildren:
        return True
      if nChildren > oChildren:
        return False
      for i in range(nChildren):
        if self.children[i] < other.children[i]:
          return True
    except AttributeError:
      return True
    return False

  def __eq__(self, other):
    return not self < other and not other < self


def _exampleCode():
  tree = TreeNode(None, 'root')
  for i in range(3):
    tree.AddChild('child %d' % i)
  print(tree)
  tree.GetChildren()[1].AddChild('grandchild')
  tree.GetChildren()[1].AddChild('grandchild2')
  tree.GetChildren()[1].AddChild('grandchild3')
  print(tree)
  tree.Pickle('save.pkl')
  print('prune')
  tree.PruneChild(tree.GetChildren()[1])
  print('done')
  print(tree)

  import copy
  tree2 = copy.deepcopy(tree)
  print('tree==tree2', tree == tree2)

  foo = [tree]
  print('tree in [tree]:', tree in foo, foo.index(tree))
  print('tree2 in [tree]:', tree2 in foo, foo.index(tree2))

  tree2.GetChildren()[1].AddChild('grandchild4')
  print('tree==tree2', tree == tree2)
  tree.Destroy()


if __name__ == '__main__':  # pragma: nocover
  _exampleCode()
