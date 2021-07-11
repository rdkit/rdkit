FAQs
======

Q1. 如何删除分子中的Atom Map ID 信息？

A1. 借助于atom的ClearProp方法可以删除指定定属性的信息

示例代码：

.. code-block:: py

  from rdkit import Chem
  smi='[CH3:1][C:2]([N:4]([C:9](=[O:10])[CH3:8])[CH2:5][CH2:6][NH2:7])=[O:3]'
  m=Chem.MolFromSmiles(smi)
  for a in m.GetAtoms():
      a.ClearProp('molAtomMapNumber')
  # [a.ClearProp('molAtomMapNumber') for a in m.GetAtoms()]
  newsmi = Chem.MolToSmiles(m)
  print(newsmi)

输出：

.. code-block:: console

  CC(=O)N(CCN)C(C)=O
  
  





           
           

    



