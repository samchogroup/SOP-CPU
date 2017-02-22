CC = g++ -O3

EFILE = ./sop.x
OBJS = ./sop.o ./random_generator.o ./global.o ./energy.o ./io.o ./params.o

sop.x: $(OBJS)
	@echo "linking ..."
	$(CC) -o $(EFILE) $(OBJS)

sop.o: ./sop.h ./random_generator.h ./global.h ./energy.h ./io.h ./params.h
	$(CC) -c ./sop.cpp -o ./sop.o

random_generator.o: ./random_generator.h
	$(CC) -c ./random_generator.cpp -o ./random_generator.o

global.o: ./global.h ./random_generator.h
		$(CC) -c ./global.cpp -o ./global.o

energy.o: ./global.h ./energy.h
		$(CC) -c ./energy.cpp -o ./energy.o

io.o: ./global.h ./io.h
		$(CC) -c ./io.cpp -o ./io.o

params.o: ./params.h ./params.h
		$(CC) -c ./params.cpp -o ./params.o

clean:
	rm -f $(OBJS) $(EFILE)
