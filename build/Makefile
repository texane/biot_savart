##
## Made by fabien le mentec <texane@gmail.com>
## 

CC	:=	g++
LD	:=	g++

SDL_CFLAGS :=	$(shell pkg-config --cflags sdl)
SDL_LFLAGS :=	$(shell pkg-config --libs sdl)

CFLAGS	:=	-Wall -O3 $(SDL_CFLAGS) -I../src
LFLAGS	:=	$(SDL_LFLAGS) -lm

RM	:=	rm -f

BIN	:=	bs

SRCS	:=	../src/main.cc			\
		../src/x.cc

OBJS	:=	$(SRCS:.cc=.o)


all	:	$(BIN)
$(BIN)	:	$(OBJS)
		$(LD) $(LFLAGS) -o $(BIN) $(OBJS)
clean	:
		$(RM) $(OBJS)
fclean	:	clean
		$(RM) $(BIN)
re	:	fclean all
%.o	:	%.cc
		$(CC) $(CFLAGS) -o $@ -c $<
