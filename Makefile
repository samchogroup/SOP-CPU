CC = g++ -O3 -g -O0

EFILE = ./sop.x
OBJS = ./bin/sop.o ./bin/random_generator.o ./bin/global.o ./bin/energy.o ./bin/io.o ./bin/params.o ./bin/neighbor_list.o ./bin/cell_list.o ./bin/barnes_hut.o

sop.x: $(OBJS)
	@echo "linking ..."
	$(CC) -o $(EFILE) $(OBJS)

bin/sop.o: ./sop.h ./random_generator.h ./global.h ./energy.h ./io.h ./params.h ./neighbor_list.h ./cell_list.h ./barnes_hut.h
	$(CC) -c ./sop.cpp -o ./bin/sop.o

bin/random_generator.o: ./random_generator.h
	$(CC) -c ./random_generator.cpp -o ./bin/random_generator.o

bin/global.o: ./global.h ./random_generator.h
		$(CC) -c ./global.cpp -o ./bin/global.o

bin/energy.o: ./global.h ./energy.h
		$(CC) -c ./energy.cpp -o ./bin/energy.o

bin/io.o: ./global.h ./io.h
		$(CC) -c ./io.cpp -o ./bin/io.o

bin/params.o: ./global.h ./params.h
		$(CC) -c ./params.cpp -o ./bin/params.o

bin/neighbor_list.o: ./global.h ./neighbor_list.h
		$(CC) -c ./neighbor_list.cpp -o ./bin/neighbor_list.o

bin/cell_list.o: ./global.h ./cell_list.h
		$(CC) -c ./cell_list.cpp -o ./bin/cell_list.o

bin/barnes_hut.o: ./global.h ./barnes_hut.h
		$(CC) -c ./barnes_hut.cpp -o ./bin/barnes_hut.o

clean:
	rm -f $(OBJS) $(EFILE)
