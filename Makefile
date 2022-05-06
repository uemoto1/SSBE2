include make.inc

TARGET = SSBE 

OBJS = \
math_constants.o \
io/salmon_file.o \
input_parameter.o \
em_field.o \
ground_state.o \
time_evolution.o \
main.o

# misc/unusedvar.o \
# common/structures.o \
# common/pack_unpack.o \
# sbe_gs.o \
# sbe_bloch_solver.o \
# test.o \


$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(FLAGS) $(LIBS)

.SUFFIXES: .f90

%.o: %.f90
	$(FC) -c -o $@ $^ $(FLAGS)

.PHONY: all clean 

all: $(TARGET)

clean:
	rm $(TARGET) $(OBJS) *.mod

