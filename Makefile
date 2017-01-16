
test: test.o voronoi.o
	clang++ $^ -o $@ -std=c++11 -g

%.o: %.cpp geometry.h debug.h
	clang++ $< -c -o $@ -std=c++11 -g

clean:
	rm -f test.o voronoi.o test
