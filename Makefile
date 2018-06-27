CC = g++ -O3

EFILE = ./sop.x
OBJS = ./sop.o ./random_generator.o ./global.o ./energy.o ./io.o ./params.o ./neighbor_list.o ./cell_list.o ./cell_array.o ./two_cells.o

sop.x: $(OBJS)
	@echo "linking ..."
	$(CC) -o $(EFILE) $(OBJS)

sop.o: ./sop.h ./random_generator.h ./global.h ./energy.h ./io.h ./params.h ./neighbor_list.h ./cell_list.o ./cell_array.o ./two_cells.o
	$(CC) -c ./sop.cpp -o ./sop.o

random_generator.o: ./random_generator.h
	$(CC) -c ./random_generator.cpp -o ./random_generator.o

global.o: ./global.h ./random_generator.h
		$(CC) -c ./global.cpp -o ./global.o

energy.o: ./global.h ./energy.h
		$(CC) -c ./energy.cpp -o ./energy.o

io.o: ./global.h ./io.h
		$(CC) -c ./io.cpp -o ./io.o

params.o: ./global.h ./params.h
		$(CC) -c ./params.cpp -o ./params.o

neighbor_list.o: ./global.h ./neighbor_list.h
		$(CC) -c ./neighbor_list.cpp -o ./neighbor_list.o

cell_list.o: ./global.h ./cell_list.h
		$(CC) -c ./cell_list.cpp -o ./cell_list.o

cell_array.o: ./global.h ./cell_array.h
		$(CC) -c ./cell_array.cpp -o ./cell_array.o

two_cells.o: ./global.h ./two_cells.h
		$(CC) -c ./two_cells.cpp -o ./two_cells.o

clean:
	rm -f $(OBJS) $(EFILE)
