sim: main.cpp
	g++ main.cpp -llapacke -o sim

clean:
	rm sim
