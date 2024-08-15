CXX=g++
CXXFLAGS=-I. -w  -Wall -std=c++11

DEPS = option_class.hpp utility_functions.hpp
OBJ = main.o exotic_options.o vanilla_options.o

main: $(OBJ)
	$(CXX) -o $@ $^

%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean
clean:
	rm -f *.o main
