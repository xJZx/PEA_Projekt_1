#include <string>

class TravellingSalesmanProblem
{
private:
	int V, E;
	int** array;
	
public:
	TravellingSalesmanProblem(int, int, int);
	TravellingSalesmanProblem(std::string src);
	~TravellingSalesmanProblem();

	void bruteForce();

	void bruteForce_test();

	void print();
};