#include <string>
#include <vector>

class TravellingSalesmanProblem
{
private:
	// ilosc wierzcholkow
	int V;
	// macierz s�siedztwa, graf pe�ny, skierowany
	std::vector<std::vector<int>> matrix;

	// macierz wag dla dynamicznego
	std::vector<std::vector<int>> dp;
	// macierz trasy dla dynamicznego
	std::vector<int> pathDynamic;
	
public:
	TravellingSalesmanProblem(int, int, int);
	TravellingSalesmanProblem(std::string src);
	~TravellingSalesmanProblem();

	void bruteForce();
	void littleAlgorithm();
	void dynamicProgramming();
	int dynamicProgrammingRecursion(int, int);

	void bruteForce_test();
	void littleAlgorithm_test();
	void dynamicProgramming_test();

	void print();
	void saveToFile();
};