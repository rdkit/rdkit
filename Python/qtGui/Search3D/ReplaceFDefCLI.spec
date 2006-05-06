import RDConfig
a = Analysis([os.path.join(HOMEPATH,'support\\_mountzlib.py'), os.path.join(HOMEPATH,'support\\useUnicode.py'), 'ReplaceFDefCLI.py'],
             #pathex=['c:\\glandrum\\RD\\trunk\\Python\\qtGui\\Search3D'],
	     excludes=['Tkinter','reportlab','_tkinter','gvib',
	     'pyPgSQL','PIL._imagingtk','_imagingtk','qt'])

pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=1,
          name='buildReplaceFDefCLI/ReplaceFDefCLI.exe',
          debug=0,
          strip=0,
          upx=0,
          console=1,
)
coll = COLLECT( exe,
               a.binaries,
               strip=0,
               upx=0,
               name='distReplaceFDefCLI')
