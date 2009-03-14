# derived from openbabel's openbabel_postinstall.py script
import os, sys
import _winreg

if len(sys.argv)==2 and sys.argv[1]=="-install":
    # Connect to the registry
    registry = _winreg.ConnectRegistry(None,_winreg.HKEY_LOCAL_MACHINE)
    # Open Environment key for writing
    environment =_winreg.OpenKey(registry,
       r"SYSTEM\CurrentControlSet\Control\Session Manager\Environment",
       0,_winreg.KEY_ALL_ACCESS)
    # Set the value of RDBASE
    rdbasedir = os.path.join(sys.prefix, "share", "rdkit")
    _winreg.SetValueEx(environment, "RDBASE", 0, _winreg.REG_EXPAND_SZ,
       rdbasedir)
    _winreg.CloseKey(environment)
    _winreg.CloseKey(registry)

    print "RDBASE is set to %s" % rdbasedir
    print
    print "You will need to reboot before the new value of RDBASE"
    print "is available, but the module should mostly work without"
    print "a reboot."
