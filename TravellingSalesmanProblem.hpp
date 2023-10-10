#include <string>

class TravellingSalesmanProblem
{
private:
	// ilosc wierzcholkow
	int V, E;
	// macierz s�siedztwa, graf pe�ny, skierowany
	int** array;
	
public:
	TravellingSalesmanProblem(int, int, int);
	TravellingSalesmanProblem(std::string src);
	~TravellingSalesmanProblem();

	void bruteForce();
	void littleAlgorithm();

	void bruteForce_test();

	void print();
};