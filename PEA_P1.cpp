// PEA_P1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "TravellingSalesmanProblem.hpp"
#include "TravellingSalesmanProblemBenchmark.hpp"
#include <iostream>

int main()
{
    while (true) {
        int mode = 0;
        std::cout << "Select the mode:\n";
        std::cout << "1. Normal\n";
        std::cout << "2. Benchmark\n";
        std::cout << "Any button - EXIT\n";
        std::cin >> mode;

        if (mode == 1) {
            int creation = 0;
            std::cout << "Select how the graph will be created:\n";
            std::cout << "1. Random (NOT WORKING YET)\n";
            std::cout << "2. From the file\n";
            std::cin >> creation;
            switch (creation) {
            case 1:
            {
                int v;
                int low, high;
                std::cout << "Insert the amount of vertices: ";
                std::cin >> v;
                std::cout << "Insert lower range: ";
                std::cin >> low;
                std::cout << "Insert upper range: ";
                std::cin >> high;
                std::cout << std::endl;

                TravellingSalesmanProblem* travellingSalesmanProblem = new TravellingSalesmanProblem(v, low, high);
                travellingSalesmanProblem->print();
                travellingSalesmanProblem->bruteForce();
                travellingSalesmanProblem->littleAlgorithm();

                delete travellingSalesmanProblem;
            }
            break;

            case 2:
            {
                std::string filename;
                std::cout << "Insert the filename: ";
                std::cin >> filename;
                std::cout << std::endl;
                TravellingSalesmanProblem* travellingSalesmanProblem = new TravellingSalesmanProblem(filename);
                travellingSalesmanProblem->print();
                travellingSalesmanProblem->bruteForce();
                travellingSalesmanProblem->littleAlgorithm();

                delete travellingSalesmanProblem;
            }
            break;

            }
        }
        else if (mode == 2) {
            int creation = 0;
            std::cout << "Select how the graph will be created:\n";
            std::cout << "1. Random (NOT WORKING YET)\n";
            std::cout << "2. From the file\n";
            std::cin >> creation;
            switch (creation) {
            case 1:
            {
                int v;
                int low, high;
                std::cout << "Insert the amount of vertices: ";
                std::cin >> v;
                std::cout << "Insert lower range: ";
                std::cin >> low;
                std::cout << "Insert upper range: ";
                std::cin >> high;
                std::cout << std::endl;

                TravellingSalesmanProblemBenchmark* travellingSalesmanProblemBenchmark = new TravellingSalesmanProblemBenchmark(v, low, high);
                delete travellingSalesmanProblemBenchmark;
            }
            break;

            case 2:
            {
                std::string filename;
                std::cout << "Insert the filename: ";
                std::cin >> filename;
                std::cout << std::endl;
                TravellingSalesmanProblemBenchmark* travellingSalesmanProblemBenchmark = new TravellingSalesmanProblemBenchmark(filename);

                delete travellingSalesmanProblemBenchmark;
            }
            break;
            }
        }
        else {
            return 0;
        }

    }
    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
