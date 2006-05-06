import RDConfig
a = Analysis([os.path.join(HOMEPATH,'support\\_mountzlib.py'), os.path.join(HOMEPATH,'support\\useUnicode.py'), 'Search3D.py'],
             pathex=['c:\\glandrum\\RD\\trunk\\Python\\qtGui\\Search3D'],
	     excludes=['Tkinter','reportlab','_tkinter','gvib',
	     'pyPgSQL','PIL._imagingtk','_imagingtk'])

dataDir=  os.path.join(RDConfig.RDDocsDir,'Programs','RDPharm3D')
dataFiles = TOC((\
('Docs/Programs/RDPharm3D/RDPharm3D-Splash.jpg', \
    os.path.join(dataDir,'RDPharm3D-Splash.jpg'),
    'DATA'),
('Docs/Programs/RDPharm3D/RDPharm3D.chm', \
    os.path.join(dataDir,'RDPharm3D.chm'),
    'DATA'),
))
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=1,
          name='buildSearch3D/RDPharm3D.exe',
          debug=0,
          strip=0,
          upx=0,
          console=1,
	  icon=os.path.join(dataDir,'Icon.ico'),
)
coll = COLLECT( exe,
               a.binaries+dataFiles,
               strip=0,
               upx=0,
               name='distSearch3D')
