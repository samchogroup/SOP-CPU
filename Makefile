CC = g++ -O3

EFILE = ./sop.x
OBJS = ./sop.o ./random_generator.o

sop.x: $(OBJS)	
	@echo "linking ..."	
	$(CC) -o $(EFILE) $(OBJS)

sop.o: ./sop.h ./random_generator.h
	$(CC) -c ./sop.cpp -o ./sop.o

random_generator.o: ./random_generator.h
	$(CC) -c ./random_generator.cpp -o ./random_generator.o

clean:
	rm -f $(OBJS) $(EFILE)
