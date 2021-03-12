CXX = g++
CXXFLAGS = -O2 -Wall -std=c++17 -Iinclude
LDFLAGS = 

EXE = matrix
TEST = test
SRCDIR = src
BINDIR = obj

OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(BINDIR)/%.o,$(wildcard $(SRCDIR)/*.cpp))

all: $(EXE)

$(EXE): $(BINDIR) $(OBJECTS) obj/main.o
	$(CXX) $(OBJECTS) obj/main.o -o $(EXE) $(LDFLAGS)

$(TEST): $(BINDIR) $(OBJECTS) obj/test.o obj/test_class.o
	$(CXX) $(OBJECTS) obj/test.o obj/test_class.o -o $(TEST) $(LDFLAGS)

obj/main.o: src/main/main.cpp
	$(CXX) $(CXXFLAGS) -c src/main/main.cpp -o obj/main.o

obj/test.o: tests/test.cpp
	$(CXX) $(CXXFLAGS) -c tests/test.cpp -o obj/test.o

obj/test_class.o: tests/test_class.cpp include/tests/test_class.hpp
	$(CXX) $(CXXFLAGS) -c tests/test_class.cpp -o obj/test_class.o

$(BINDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -MMD -o $@ $<

include $(wildcard $(BINDIR)/*.d)

$(BINDIR):
	mkdir -p $(BINDIR)

clean:
	rm -rf $(BINDIR) $(EXE) $(TEST)
.PHONY: clean all

