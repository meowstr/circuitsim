all: tran sim

tran: tran.cpp
	g++ tran.cpp -o tran

sim: main.cpp
	g++ main.cpp -llapacke -o sim

clean:
	rm sim
	rm tran
