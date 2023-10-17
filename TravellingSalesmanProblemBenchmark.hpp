#include <string>

class TravellingSalesmanProblemBenchmark
{
	private:
		void test_bruteForce(std::string src);
		void test_littleAlgorithm(std::string src);

	public:
		TravellingSalesmanProblemBenchmark(std::string src, std::string type);
};