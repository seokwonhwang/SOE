CPP = g++

EXEC = ejection

#PERROR = -DPERROR

PRINT = -DPRINT

CFLAGS = -O3 -Wall $(PRINT) $(PERROR)

HOMEDIR = /home/wwalker/UPS/Materials/

INCLUDE = -I$(HOMEDIR) \

LIBS =  -lm

SRCS =  ejection_test.cpp Ejection.cpp

OBJS = $(SRCS:.cpp=.o)

all : $(EVAL) $(LOAD) ${EXEC}

$(EXEC): $(OBJEVAL) $(OBJLOAD) $(OBJS)
	$(CPP) -o $(EXEC) $(OBJS) $(OBJEVAL) $(OBJLOAD) $(LIBS) 

.cpp.o:
	$(CPP) $(CFLAGS) $(INCLUDE) -c $< 
clean:
	rm -f $(OBJS) $(EXEC)
