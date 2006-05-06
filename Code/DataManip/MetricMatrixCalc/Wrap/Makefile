
$(RDBASE)/Python/DataManip/Metric/rdMetricMatrixCalc.so: rdMetricMatrixCalc.cpp
	python setup.py install --install-lib=$(RDBASE)/Python

clean:
	if test -e $(RDBASE)/Python/DataManip/Metric/rdMetricMatrixCalc.so; then rm $(RDBASE)/Python/DataManip/Metric/rdMetricMatrixCalc.so; fi
	if test -e build; then rm -rf build; fi
