"""Patch rdkit/__init__.py to auto-load .so.wasm files on emscripten.

This makes 'micropip.install()' + 'import rdkit' just work:
1. .so files are renamed to .so.wasm so micropip doesn't auto-load them
2. __init__.py loads librdkit_core.so.wasm with {global: true} on first import
3. A custom import finder loads wrapper .so.wasm via ExtensionFileLoader
"""
import sys

init = open("/out/rdkit/__init__.py").read()

loader = '''import sys as _sys

if _sys.platform == 'emscripten':
    import os as _os
    import importlib.abc
    import importlib.machinery

    from pyodide_js._module import loadDynamicLibrary as _ldl
    import js as _js
    _core = _os.path.join(_os.path.dirname(__file__), 'librdkit_core.so.wasm')
    _ldl(_core, _js.JSON.parse('{"global": true}'))
    del _ldl, _js, _core

    class _RDKitExtensionFinder(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path, target=None):
            parts = fullname.split('.')
            if parts[0] != 'rdkit':
                return None
            modname = parts[-1]
            if path:
                for d in path:
                    candidate = _os.path.join(d, modname + '.so.wasm')
                    if _os.path.exists(candidate):
                        loader = importlib.machinery.ExtensionFileLoader(
                            fullname, candidate
                        )
                        return importlib.util.spec_from_file_location(
                            fullname, candidate, loader=loader,
                        )
            return None

    import importlib.util
    _sys.meta_path.insert(0, _RDKitExtensionFinder())

'''

open("/out/rdkit/__init__.py", "w").write(loader + init)
print("Patched __init__.py")
