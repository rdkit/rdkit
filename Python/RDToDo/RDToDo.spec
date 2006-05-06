a = Analysis(['c:/installer\\support\\_mountzlib.py', 'c:\\glandrum\\RD\\trunk\\Python\\RDToDo\\RDToDo.py'],
             pathex=[],
				excludes=['PIL','Tkinter','win32com'])

pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=1,
          name='buildRDToDo/RDToDo.exe',
          debug=0,
          strip=0,
          console=1 )
coll = COLLECT( exe,
               a.binaries,
               strip=0,
               name='distRDToDo')
