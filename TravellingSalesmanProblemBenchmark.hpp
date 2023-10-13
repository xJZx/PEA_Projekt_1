#include <string>

class TravellingSalesmanProblemBenchmark
{
	private:
		void test_bruteForce(int v, int low, int high);
		void test_littleAlgorithm(int v, int low, int high);

		void test_bruteForce(std::string src);
		void test_littleAlgorithm(std::string src);

	public:
		TravellingSalesmanProblemBenchmark(int v, int low, int high, std::string type);
		TravellingSalesmanProblemBenchmark(std::string src, std::string type);
};