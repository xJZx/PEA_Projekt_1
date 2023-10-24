#include <string>
#include <vector>

class TravellingSalesmanProblem
{
private:
	// ilosc wierzcholkow
	int V;
	// macierz s¹siedztwa, graf pe³ny, skierowany
	std::vector<std::vector<int>> matrix;
	
public:
	TravellingSalesmanProblem(int, int, int);
	TravellingSalesmanProblem(std::string src);
	~TravellingSalesmanProblem();

	void bruteForce();
	void littleAlgorithm();
	void dynamicProgramming();
	int dynamicProgrammingRecursion(int, int, std::vector<std::vector<int>>, std::vector<int>);

	void bruteForce_test();
	void littleAlgorithm_test();
	void dynamicProgramming_test();

	void print();
	void saveToFile();
};