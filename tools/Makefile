.PHONY: all clean

CC := h5cc
LDFLAGS := -L ../external_libraries/silo/lib
LDLIBS := -lsiloh5 -lm
CPPFLAGS := -I ../external_libraries/silo/include -Wall

all: silo_to_raw

clean:
	$(RM) -f silo_to_raw silo_to_raw.o
