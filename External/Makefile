
all:
	$(MAKE) --directory=$(RDBASE)/Code/RDGeneral install
	$(MAKE) --directory=$(RDBASE)/Code/RDBoost
	$(MAKE) --directory=Lapack++ libs
	$(MAKE) --directory=svdlibc libsvd.a
	cd svdlibc/Wrap && python setup.py install --install-lib=$(RDBASE)/Python
	$(MAKE) --directory=svdpackc sis2-rd && cp svdpackc/sis2-rd $(RDBASE)/bin
	$(MAKE) --directory=vflib-2.0
	$(MAKE) --directory=cmim-1.0 libfastentropy.a
