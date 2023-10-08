#include "TravellingSalesmanProblemBenchmark.hpp"
#include "TravellingSalesmanProblem.hpp"
#include <time.h>
#include <windows.h>
#include <iostream>

long long int read_QPC()
{
    LARGE_INTEGER count;
    QueryPerformanceCounter((LARGE_INTEGER*)&count);
    return ((long long int) count.QuadPart);
}

TravellingSalesmanProblemBenchmark::TravellingSalesmanProblemBenchmark(std::string src)
{
    test_bruteForce(src);
}

void TravellingSalesmanProblemBenchmark::test_bruteForce(std::string src) {

    int results[5];

    for (int i = 0; i < 5; i++)
    {
        TravellingSalesmanProblem* travellingSalesmanProblem = new TravellingSalesmanProblem(src);

        long long int frequency, start, end, elapsed;
        QueryPerformanceFrequency((LARGE_INTEGER*)&frequency);

        start = read_QPC();
        travellingSalesmanProblem->bruteForce_test();
        end = read_QPC();

        elapsed = end - start;
        // przez 1000, aby byly milisekundy
        results[i] = 1000.0 * elapsed / frequency;
    }
    int result = 0;
    for (int i = 0; i < 5; i++)
    {
        result += results[i];
    }

    std::cout << "TSP Brute Force algorithm average time: " << result / 5 << "[ms]" << std::endl;
}