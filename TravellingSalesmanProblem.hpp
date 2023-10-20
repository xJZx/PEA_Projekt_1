#include <string>
#include <vector>

class TravellingSalesmanProblem
{
private:
	// ilosc wierzcholkow
	int V;
	// macierz s�siedztwa, graf pe�ny, skierowany
	std::vector<std::vector<int>> array;
	
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