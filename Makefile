SRCS := dbl_pnd.cpp
OBJ := $(patsubst %.cpp,%.o,$(SRCS))
BIN := $(patsubst %.cpp,%.bin,$(SRCS))
CXX := g++
CXXFLAGS := -Wall -fexceptions -pedantic -Wextra -fopenmp -Ofast
LD := g++
LDFLAGS := -Ofast -flto -fopenmp

$(BIN): $(OBJ)
	$(LD) $(LDFLAGS) $^ -o $@

clean:
	-$(RM) *.o *.bin
