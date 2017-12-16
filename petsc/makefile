CC = mpif90
CFLAGS = -O3 -Wall -I$(PETSC_DIR)/include
LDFLAGS = -lpetsc -L$(PETSC_DIR)/lib -lm
OS := $(shell uname)
ifeq ($(OS),Linux)
	LDFLAGS += -Wl,-rpath=$(PETSC_DIR)/lib
endif

TARGET = ex11f90

all: $(TARGET)

$(TARGET): $(TARGET).F90
	$(CC) $(CFLAGS) -o $(TARGET) $< $(LDFLAGS)

clean:
	rm -f *.o $(TARGET)
