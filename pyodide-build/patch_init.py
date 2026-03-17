"""Patch rdkit/__init__.py to auto-load librdkit_core.so on emscripten."""

init = open('/out/rdkit/__init__.py').read()

loader = '''import os as _os, sys as _sys
if _sys.platform == 'emscripten':
    from pyodide_js._module import loadDynamicLibrary as _ldl
    import js as _js
    _core = _os.path.join(_os.path.dirname(__file__), 'librdkit_core.so')
    _r = _ldl(_core, _js.JSON.parse('{"global": true}'))
    del _ldl, _js, _core, _r

'''

open('/out/rdkit/__init__.py', 'w').write(loader + init)
print('Patched __init__.py with core loader')
