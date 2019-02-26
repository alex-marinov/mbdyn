CFLAGS = -c -Wall -O2 -g

default: Post

lex.Mbdynpost.cc: Mbdynpost.flex Post.hh
	flex++ Mbdynpost.flex

lex.Mbdynpost.o: lex.Mbdynpost.cc
	g++ $(CFLAGS) lex.Mbdynpost.cc

Post.o: Post.C Post.hh
	g++ $(CFLAGS) Post.C

Post: lex.Mbdynpost.o Post.o
	g++ -o Post -g -Wall lex.Mbdynpost.o Post.o

clean:
	rm -f *.o Post lex.Mbdynpost.cc
