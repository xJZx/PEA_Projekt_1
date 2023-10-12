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
        // po ostatnim wierzcholku dodajemy kolejny, potem potrzebne do permutacji
        trip.insert(trip.end(), i);
    }

    // okreslenie powrotu do punktu wyjscia dla domyslnej sciezki
    trip.insert(trip.end(), source);

    // next_permutation z biblioteki <algorithm>, pozwala rearan¿owaæ elementy tab. wekt. tak aby uzyskaæ wszystkie mo¿liwe jej permutacje
    // rearan¿ujemy tylko wierzcho³ki odwiedzane, nie wyjœciowy
    do {
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
    } while (std::next_permutation(trip.begin() + 1, trip.end() - 1));

    // V + 1, aby pokazaæ ¿e komiwoja¿er dotar³ do domku
    std::cout << "Best trip: ";
    for (int i = 0; i < V + 1; i++) {
        std::cout << bestTrip[i] << " ";
    }

    std::cout << std::endl;
    std::cout << "Best value: " << bestValue << std::endl;
}

void TravellingSalesmanProblem::bruteForce_test() {
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
}

void TravellingSalesmanProblem::littleAlgorithm() {
    // tablica wspolczynnikow standaryzacji
    std::vector<int> aFactor;
    std::vector<int> bFactor;

    int min = INT_MAX;

    int lowerBound = 0;

    for (int N = V; N > 0; N--) {
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                // szukamy minimum dla RZÊDU
                if (array[row][column] < min && row != column && array[row][column] != -1) {
                    min = array[row][column];
                }
            }
            // dodajemy minimum do tablicy wspó³cz. a
            aFactor.emplace_back(min);
            std::cout << "a" << row << " " << min << std::endl;
            // ponowne ustalenie minimum
            min = INT_MAX;
        }

        // 1. krok do stworzenia C' - odjêcie od C wspó³czynnika a
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                // Cij - ai
                if (row != column && array[row][column] != -1) {
                    array[row][column] -= aFactor.at(row);
                }
            }
        }

        //print();

        for (int column = 0; column < V; column++) {
            for (int row = 0; row < V; row++) {
                // szukamy minimum dla KOLUMNY
                if (array[row][column] < min && row != column && array[row][column] != -1) {
                    min = array[row][column];
                }
            }
            // dodajemy minimum do tablicy wspó³cz. a
            bFactor.emplace_back(min);
            std::cout << "b" << column << " " << min << std::endl;
            // ponowne ustalenie minimum
            min = INT_MAX;
        }

        // 2. krok do stworzenia C' - odjêcie od C wspó³czynnika b
        for (int column = 0; column < V; column++) {
            for (int row = 0; row < V; row++) {
                // Cij - bi
                if (column != row && array[column][row] != -1) {
                    array[column][row] -= bFactor.at(row);
                }
            }
        }

        //print();

        // obliczenie dolnego oszacowania dla wszystkich rozwi¹zañ
        for (int i = 0; i < V; i++) {
            if (aFactor[i] == INT_MAX) {
                lowerBound += 0;
            }
            else {
                lowerBound += aFactor[i];
            }
            if (bFactor[i] == INT_MAX) {
                lowerBound += 0;
            }
            else {
                lowerBound += bFactor[i];
            }
        }

        std::cout << "Lower bound on all tours: " << lowerBound << std::endl;

        // wyczyszczenie tablic wektorowych
        aFactor.clear();
        bFactor.clear();

        // stworzenie macierzy kosztów rezygnacji dla tras "zerowych"
        std::vector<std::vector<int>> resignationArray(V, std::vector<int>(V, -1));
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                if (array[row][column] == 0 && row != column && array[row][column] != -1) {
                    int minRow = INT_MAX;
                    int minColumn = INT_MAX;
                    for (int i = 0; i < V; i++) {
                        if (array[row][i] < minRow && row != i && i != column && array[row][i] != -1) {
                            minRow = array[row][i];
                        }
                    }
                    for (int j = 0; j < V; j++) {
                        if (array[j][column] < minColumn && j != column && j != row && array[j][column] != -1) {
                            minColumn = array[j][column];
                        }
                    }
                    if (minRow != INT_MAX && minColumn != INT_MAX) {
                        resignationArray[row][column] = minRow + minColumn;
                    }
                }
            }
        }

        // sprawdzenie poprawnosci macierzy rezygnacji
        std::cout << "Resignation Matrix: " << std::endl;
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                std::cout << resignationArray[i][j] << " ";
            }
            std::cout << std::endl;
        }

        // znalezienie d(kl) dla macierzy rezygnacji
        int d_kl = 0;
        // od razu zapisujemy row i column które bêd¹ wykreœlane z macierzy
        int k, l;
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                if (resignationArray[row][column] >= d_kl) {
                    d_kl = resignationArray[row][column];
                    k = row;
                    l = column;
                }
            }
        }
        std::cout << "d_kl = " << d_kl << std::endl;

        // sprawdzenie warunku koñcowego d_kl == 0
        if (d_kl == 0) {
            return;
        }

        // obliczenie dolnej granicy dla K2
        int lowerBound_K2 = lowerBound + d_kl;

        // tworzymy macierz zredukowan¹ C1
        // najpierw usuwamy rz¹d
        for (int i = 0; i < V; i++) {
            array[k][i] = -1;
        }
        // potem kolumnê
        for (int j = 0; j < V; j++) {
            array[j][l] = -1;
        }
        // blokujemy podcykl tej samej œcie¿ki
        array[l][k] = -1;

        print();
    }

}

void TravellingSalesmanProblem::littleAlgorithm_test() {
    // tablica wspolczynnikow standaryzacji
    std::vector<int> aFactor;
    std::vector<int> bFactor;

    int min = INT_MAX;

    int lowerBound = 0;

    for (int N = V; N > 0; N--) {
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                // szukamy minimum dla RZÊDU
                if (array[row][column] < min && row != column && array[row][column] != -1) {
                    min = array[row][column];
                }
            }
            // dodajemy minimum do tablicy wspó³cz. a
            aFactor.emplace_back(min);
            // ponowne ustalenie minimum
            min = INT_MAX;
        }

        // 1. krok do stworzenia C' - odjêcie od C wspó³czynnika a
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                // Cij - ai
                if (row != column && array[row][column] != -1) {
                    array[row][column] -= aFactor.at(row);
                }
            }
        }

        for (int column = 0; column < V; column++) {
            for (int row = 0; row < V; row++) {
                // szukamy minimum dla KOLUMNY
                if (array[row][column] < min && row != column && array[row][column] != -1) {
                    min = array[row][column];
                }
            }
            // dodajemy minimum do tablicy wspó³cz. a
            bFactor.emplace_back(min);
            // ponowne ustalenie minimum
            min = INT_MAX;
        }

        // 2. krok do stworzenia C' - odjêcie od C wspó³czynnika b
        for (int column = 0; column < V; column++) {
            for (int row = 0; row < V; row++) {
                // Cij - bi
                if (column != row && array[column][row] != -1) {
                    array[column][row] -= bFactor.at(row);
                }
            }
        }

        // obliczenie dolnego oszacowania dla wszystkich rozwi¹zañ
        for (int i = 0; i < V; i++) {
            if (aFactor[i] == INT_MAX) {
                lowerBound += 0;
            }
            else {
                lowerBound += aFactor[i];
            }
            if (bFactor[i] == INT_MAX) {
                lowerBound += 0;
            }
            else {
                lowerBound += bFactor[i];
            }
        }

        // wyczyszczenie tablic wektorowych
        aFactor.clear();
        bFactor.clear();

        // stworzenie macierzy kosztów rezygnacji dla tras "zerowych"
        std::vector<std::vector<int>> resignationArray(V, std::vector<int>(V, -1));
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                if (array[row][column] == 0 && row != column && array[row][column] != -1) {
                    int minRow = INT_MAX;
                    int minColumn = INT_MAX;
                    for (int i = 0; i < V; i++) {
                        if (array[row][i] < minRow && row != i && i != column && array[row][i] != -1) {
                            minRow = array[row][i];
                        }
                    }
                    for (int j = 0; j < V; j++) {
                        if (array[j][column] < minColumn && j != column && j != row && array[j][column] != -1) {
                            minColumn = array[j][column];
                        }
                    }
                    if (minRow != INT_MAX && minColumn != INT_MAX) {
                        resignationArray[row][column] = minRow + minColumn;
                    }
                }
            }
        }

        // znalezienie d(kl) dla macierzy rezygnacji
        int d_kl = 0;
        // od razu zapisujemy row i column które bêd¹ wykreœlane z macierzy
        int k, l;
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                if (resignationArray[row][column] >= d_kl) {
                    d_kl = resignationArray[row][column];
                    k = row;
                    l = column;
                }
            }
        }

        // sprawdzenie warunku koñcowego d_kl == 0
        if (d_kl == 0) {
            return;
        }

        // obliczenie dolnej granicy dla K2
        int lowerBound_K2 = lowerBound + d_kl;

        // tworzymy macierz zredukowan¹ C1
        // najpierw usuwamy rz¹d
        for (int i = 0; i < V; i++) {
            array[k][i] = -1;
        }
        // potem kolumnê
        for (int j = 0; j < V; j++) {
            array[j][l] = -1;
        }
        // blokujemy podcykl tej samej œcie¿ki
        array[l][k] = -1;
    }
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