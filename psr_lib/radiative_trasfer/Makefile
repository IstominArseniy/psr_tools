CC=g++ 
CFLAGS=-std=c++17 -I /usr/local/include/eigen #path to Eigen3 library
LIBDIR=cpp_lib
BINDINGDIR=cpp_lib
LIBSOURCES=$(LIBDIR)/CppInterfaceClass.cpp $(LIBDIR)/SolverClass.cpp $(LIBDIR)/PSRClass.cpp $(LIBDIR)/FixedHeightModelClass.cpp $(LIBDIR)/constants.cpp $(LIBDIR)/functions.cpp $(LIBDIR)/integrator.cpp $(LIBDIR)/read_write.cpp
BINDINGSOURCES=$(BINDINGDIR)/interface_class_binding.cpp 


all: python_binding

python_binding: $(BINDINGSOURCES) $(LIBSOURCES)
	$(CC) $(CFLAGS) -shared -fPIC $$(python3 -m pybind11 --includes) $(BINDINGSOURCES) $(LIBSOURCES) -o ./cpp_interface$$(python3-config --extension-suffix)

clean:
	rm -rf $(OBJDIR)/*.o
