Script downloaded from https://github.com/x2on/expat-ios

Patch from http://stackoverflow.com/questions/10453946/how-to-compile-expat-with-ios-sdk-5-1

To build with the current directory layout in expatpp, 

(gave up, currently copying everything into my projects)
- copy in the contents of expat and expatpp
- add the expat_config.h defined by running configure
- define HAVE_EXPAT_CONFIG_H
- change the #includes of expat_config.h from <> to "" so the local one found
