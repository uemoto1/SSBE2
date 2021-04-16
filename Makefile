include Makefile.inc

TARGET = SSBE 

OBJS = \
misc/unusedvar.o \
math/salmon_math.o \
io/salmon_file.o \
common/structures.o \
common/pack_unpack.o \
sbe_gs.o \
sbe_bloch_solver.o \
pulse.o \
test.o \
input_parameter.o \
main.o

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(FLAGS) $(LIBS)

.SUFFIXES: .f90

%.o: %.f90
	$(FC) -c -o $@ $^ $(FLAGS)

.PHONY: all clean 

all: $(TARGET)

clean:
	rm $(TARGET) $(OBJS) *.mod

input_parameter.f90: input_parameter.py
	python input_parameter.py > $@
