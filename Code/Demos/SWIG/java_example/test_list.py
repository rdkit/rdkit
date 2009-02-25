import os
classP = os.environ.get("CLASSPATH","").split(":")
if "." not in classP:
    classP.append(".")
classP.append("RDKit_Wrapper.jar")
classP = ":".join(classP)
tests=[
  ("java","-cp %s WrapperTests"%classP,{}),
  ]



longTests=[

  ]

if __name__=='__main__':
  import sys
  from rdkit import TestRunner
  failed,tests = TestRunner.RunScript('test_list.py',0,1)
  sys.exit(len(failed))
