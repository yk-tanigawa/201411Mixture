CC = g++
LD = g++
CFLAGS =  -O2 -Warning
LDFLAGS =  -lpthread
SRCS := $(wildcard *.cpp) # wildcard
OBJS = $(SRCS:.cpp=.o)
DEPS = $(SRCS:.cpp=.dep)
EXEC = $(SRCS:.cpp=)
RM = rm -f


all: EM

EM: EM.o
	$(LD)  -o $@ $^ $(LDFLAGS)

plot:
	$(RM) output/results.gif
	R --slave --vanilla < plot.R

clean:
	$(RM) $(OBJS) $(EXEC) *~

check-syntax:
	$(CC) -Wall -Wextra -pedantic -fsyntax-only $(CHK_SOURCES)

.PHONY: clean 
	all clean
