%.d: %.cpp
	@echo Making dependencies for file $<; \
	$(CXX) -M $(CXXFLAGS) $< > $@.$$$$;                      \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;     \
	rm -f $@.$$$$
OBJS=$(subst .cpp,.o,$(SOURCES))
DEPENDS=$(subst .cpp,.d,$(SOURCES))
