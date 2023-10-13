#include <string>

class TravellingSalesmanProblem
{
private:
	// ilosc wierzcholkow
	int V;
	// macierz s¹siedztwa, graf pe³ny, skierowany
	int** array;
	
public:
	TravellingSalesmanProblem(int, int, int);
	TravellingSalesmanProblem(std::string src);
	~TravellingSalesmanProblem();

	void bruteForce();
	void littleAlgorithm();

	void bruteForce_test();
	void littleAlgorithm_test();

	void print();
	void saveToFile();
};