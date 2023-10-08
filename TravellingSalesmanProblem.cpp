#include "TravellingSalesmanProblem.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

TravellingSalesmanProblem::TravellingSalesmanProblem(int v, int low, int high)
{
    this->V = v;
    // skierowany graf, w ktorym sciezki do obu tych samych wierzcholkow mog¹ byæ ró¿ne, w zale¿noœci z którego idziemy
    // czyli niesymetryczna macierz s¹siedztwa z wagami, dlatego taki wzór
    // ZMIENIÆ NA GRAF PE£NY!!!
    this->E = (int)(V * (V - 1));


    array = new int* [V];

    for (int i = 0; i < V; i++)
        array[i] = new int[V];

    // wypelnienie macierzy maksymalnymi wartosciami
    for (int i = 0; i < V; i++)
        for (int j = 0; j < V; j++)
            array[i][j] = INT_MAX;

    std::vector<int> visited;
    int prev = 0, j = 0, next, weight;

    visited.push_back(prev);

    // wypelnienie macierzy losowymi wartosciami
    for (int i = 0; i < V - 1; i++)
    {
        next = rand() % V;
        while (j < visited.size())
        {
            if (next == visited[j])
            {
                next = (next + 1) % V;
                j = 0;
            }
            else
                j++;
        }
        j = 0;
        weight = rand() % high + low;
        array[prev][next] = weight;
        prev = next;
        visited.push_back(prev);
    }
    for (int i = V - 1; i < E; i++)
    {
        prev = rand() % V;
        next = rand() % V;
        while (prev == next)
            next = rand() % V;

        weight = rand() % high + low;
        array[prev][next] = weight;
    }


}

TravellingSalesmanProblem::TravellingSalesmanProblem(std::string src)
{
    std::fstream file;

    file.open(src);

    if (file.good())
    {
        file >> V;

        array = new int* [V];

        for (int i = 0; i < V; i++)
            array[i] = new int[V];

        for (int i = 0; i < V; i++)
            for (int j = 0; j < V; j++)
                array[i][j] = INT_MAX;

        int weight;
        for (int first = 0; first < V; first++)
        {
            for (int second = 0; second < V; second++)
            {
                file >> weight;
                array[first][second] = weight;
                // bo graf asymetryczny skierowany
                //array[second][first] = weight;
            }      
        }
        file.close();
    }
    else
    {
        array = NULL;
        std::cout << "Error occurred!\n";
    }
}

TravellingSalesmanProblem::~TravellingSalesmanProblem() {
}

void TravellingSalesmanProblem::bruteForce() {
    // przyjmijmy, ¿e zaczynamy od wierzcho³ka 0
    int source = 0;
    int bestValue = INT_MAX;
    int currentValue = 0;
    // najlepsza sciezka zapamiêtana w tablicy wektorowej
    std::vector<int> bestTrip;
    // obecna sciezka zapamietana w t. wekt.
    std::vector<int> trip;

    // okreslenie domyslnej sciezki (wyjsciowej)
    for (int i = 0; i < V; i++) {
        trip.insert(trip.end(), i);
    }

    // okreslenie powrotu do punktu wyjscia dla domyslnej sciezki
    trip.insert(trip.end(), 0);

    // next_permutation z biblioteki <algorithm>, pozwala rearan¿owaæ elementy tab. wekt. tak aby uzyskaæ wszystkie mo¿liwe jej permutacje
    while (std::next_permutation(trip.begin() + 1, trip.end() - 1)) {
        currentValue = 0;
        // dodawanie wag krawêdzi grafu do ca³kowitego kosztu obecnej œcie¿ki
        for (int i = 0; i < V; i++) {
            currentValue += array[trip[i]][trip[i + 1]];
        }

        // jeœli obecna permutacja okaza³a siê korzystniejsza, to jest zapisywana jako owa
        if (currentValue < bestValue) {
            bestValue = currentValue;
            bestTrip = trip;
        }
    }

    // V + 1, aby pokazaæ ¿e komiwoja¿er dotar³ do domku
    std::cout << "Best trip: ";
    for (int i = 0; i < V + 1; i++) {
        std::cout << bestTrip[i] << " ";
    }

    std::cout << std::endl;
    std::cout << "Best value: " << bestValue << std::endl;
}


void TravellingSalesmanProblem::print()
{
    std::cout << "Graph stored in the matrix: \n";
    for (int first = 0; first < V; first++)
    {
        for (int second = 0; second < V; second++)
        {
            std::cout << std::setw(4) << array[first][second];
        }

        std::cout << "\n";
    }
    std::cout << std::endl;
}