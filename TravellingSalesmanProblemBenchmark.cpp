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

TravellingSalesmanProblemBenchmark::TravellingSalesmanProblemBenchmark(int v, int low, int high, std::string type)
{
    if (type == "bruteForce") {
        test_bruteForce(v, low, high);
    }
    else if (type == "little") {
        test_littleAlgorithm(v, low, high);
    }
}

TravellingSalesmanProblemBenchmark::TravellingSalesmanProblemBenchmark(std::string src, std::string type)
{
    if (type == "bruteForce") {
        test_bruteForce(src);
    }
    else if (type == "little") {
        test_littleAlgorithm(src);
    }
}

// iloœæ powtórzeñ
const int N = 1;

void TravellingSalesmanProblemBenchmark::test_bruteForce(int v, int low, int high) {

    int results[N];

    for (int i = 0; i < N; i++)
    {
        TravellingSalesmanProblem* travellingSalesmanProblem = new TravellingSalesmanProblem(v, low, high);

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
    for (int i = 0; i < N; i++)
    {
        result += results[i];
    }

    std::cout << "TSP Brute Force algorithm average time: " << result / N << "[ms]" << std::endl;
}

void TravellingSalesmanProblemBenchmark::test_bruteForce(std::string src) {
    
    int results[N];

    for (int i = 0; i < N; i++)
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
    for (int i = 0; i < N; i++)
    {
        result += results[i];
    }

    std::cout << "TSP Brute Force algorithm average time: " << result / N << "[ms]" << std::endl;
}

void TravellingSalesmanProblemBenchmark::test_littleAlgorithm(int v, int low, int high) {

    int results[N];

    for (int i = 0; i < N; i++)
    {
        TravellingSalesmanProblem* travellingSalesmanProblem = new TravellingSalesmanProblem(v, low, high);

        long long int frequency, start, end, elapsed;
        QueryPerformanceFrequency((LARGE_INTEGER*)&frequency);

        start = read_QPC();
        travellingSalesmanProblem->littleAlgorithm_test();
        end = read_QPC();

        elapsed = end - start;
        // przez 1000, aby byly milisekundy
        results[i] = 1000000.0 * elapsed / frequency;
    }
    int result = 0;
    for (int i = 0; i < N; i++)
    {
        result += results[i];
    }

    std::cout << "TSP Little's Algorithm average time: " << result / N << "[us]" << std::endl;
}

void TravellingSalesmanProblemBenchmark::test_littleAlgorithm(std::string src) {

    int results[N];

    for (int i = 0; i < N; i++)
    {
        TravellingSalesmanProblem* travellingSalesmanProblem = new TravellingSalesmanProblem(src);

        long long int frequency, start, end, elapsed;
        QueryPerformanceFrequency((LARGE_INTEGER*)&frequency);

        start = read_QPC();
        travellingSalesmanProblem->littleAlgorithm_test();
        end = read_QPC();

        elapsed = end - start;
        // przez 1000, aby byly milisekundy
        results[i] = 1000000.0 * elapsed / frequency;
    }
    int result = 0;
    for (int i = 0; i < N; i++)
    {
        result += results[i];
    }

    std::cout << "TSP Little's Algorithm average time: " << result / N << "[us]" << std::endl;
}