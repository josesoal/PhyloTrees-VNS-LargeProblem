CC = gcc
CFLAGS = -g -Wall # flags when *compiling*
LFLAGS = -g -Wall # flags when *linking*
LIBS = -lm # math library
SOURCES = main.c vns.c tree.c iterate_tree.c median_solvers.c queue.c convert.c cheaptsp.c bbtsp.c inversion_median_alberto.c invdist.c dcjdist.c uf.c binencode.c random.c measure_time.c auxiliary.c condense.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = main

all: $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $@ $(LIBS)

%.o:%.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
	
	
