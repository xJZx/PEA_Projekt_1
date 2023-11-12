// PEA_P1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "TravellingSalesmanProblem.hpp"
#include "TravellingSalesmanProblemBenchmark.hpp"
#include <iostream>

int main()
{
    std::string lastFilename;

    while (true) {
        int creation = 0;
        std::cout << "Select how the graph will be created:\n";
        std::cout << "Any key - EXIT\n";
        std::cout << "1. Random\n";
        std::cout << "2. From the file\n";
        std::cout << "3. Show last\n";
        std::cout << "4. Brute Force\n";
        std::cout << "5. Little's Algorithm\n";
        std::cout << "6. Dynamic Programming\n";
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
            travellingSalesmanProblem->saveToFile();
            lastFilename = "generated.txt";

            delete travellingSalesmanProblem;
        }
        break;

        case 2:
        {
            std::string filename;
            std::cout << "Insert the filename: ";
            std::cin >> filename;
            lastFilename = filename;
            std::cout << std::endl;
        }
        break;

        case 3:
        {
            TravellingSalesmanProblem* travellingSalesmanProblem = new TravellingSalesmanProblem(lastFilename);
            travellingSalesmanProblem->print();

            delete travellingSalesmanProblem;
        }
        break;

        case 4:
        {
            TravellingSalesmanProblem* travellingSalesmanProblem = new TravellingSalesmanProblem(lastFilename);
            //travellingSalesmanProblem->bruteForce();
            TravellingSalesmanProblemBenchmark* travellingSalesmanProblemBenchmark = new TravellingSalesmanProblemBenchmark(lastFilename, "bruteForce");

            delete travellingSalesmanProblem;
            delete travellingSalesmanProblemBenchmark;
        }
        break;

        case 5:
        {
            TravellingSalesmanProblem* travellingSalesmanProblem = new TravellingSalesmanProblem(lastFilename);
            //travellingSalesmanProblem->littleAlgorithm();
            TravellingSalesmanProblemBenchmark* travellingSalesmanProblemBenchmark = new TravellingSalesmanProblemBenchmark(lastFilename, "little");

            delete travellingSalesmanProblem;
            delete travellingSalesmanProblemBenchmark;
        }
        break;

        case 6:
        {
            TravellingSalesmanProblem* travellingSalesmanProblem = new TravellingSalesmanProblem(lastFilename);
            //travellingSalesmanProblem->dynamicProgramming();
            TravellingSalesmanProblemBenchmark* travellingSalesmanProblemBenchmark = new TravellingSalesmanProblemBenchmark(lastFilename, "dynamic");

            delete travellingSalesmanProblem;
            delete travellingSalesmanProblemBenchmark;
        }
        break;

        default:
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
