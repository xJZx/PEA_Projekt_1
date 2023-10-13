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
    // skierowany graf, w ktorym sciezki do obu tych samych wierzcholkow mog� by� r�ne, w zale�no�ci z kt�rego idziemy
    // czyli niesymetryczna macierz s�siedztwa z wagami, dlatego taki wz�r
    // ZMIENI� NA GRAF PE�NY!!!
    
    array = new int* [V];

    for (int i = 0; i < V; i++)
        array[i] = new int[V];

    // wypelnienie macierzy maksymalnymi wartosciami
    for (int i = 0; i < V; i++)
        for (int j = 0; j < V; j++)
            array[i][j] = -1;

    int weight;
    for (int first = 0; first < V; first++)
    {
        for (int second = 0; second < V; second++)
        {
            if (first != second) {
                weight = (rand() % (high - low + 1)) + low;
                array[first][second] = weight;
            }
            // bo graf asymetryczny skierowany
            //array[second][first] = weight;
        }
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
    if (array == NULL) {
        std::cout << "Array is empty!\n";
        return;
    }
    // przyjmijmy, �e zaczynamy od wierzcho�ka 0
    int source = 0;
    int bestValue = INT_MAX;
    int currentValue = 0;
    // najlepsza sciezka zapami�tana w tablicy wektorowej
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

    // next_permutation z biblioteki <algorithm>, pozwala rearan�owa� elementy tab. wekt. tak aby uzyska� wszystkie mo�liwe jej permutacje
    // rearan�ujemy tylko wierzcho�ki odwiedzane, nie wyj�ciowy
    do {
        currentValue = 0;
        // dodawanie wag kraw�dzi grafu do ca�kowitego kosztu obecnej �cie�ki
        for (int i = 0; i < V; i++) {
            currentValue += array[trip[i]][trip[i + 1]];
        }

        // je�li obecna permutacja okaza�a si� korzystniejsza, to jest zapisywana jako owa
        if (currentValue < bestValue) {
            bestValue = currentValue;
            bestTrip = trip;
        }
    } while (std::next_permutation(trip.begin() + 1, trip.end() - 1));

    // V + 1, aby pokaza� �e komiwoja�er dotar� do domku
    std::cout << "Best trip: ";
    for (int i = 0; i < V + 1; i++) {
        std::cout << bestTrip[i] << " ";
    }

    std::cout << std::endl;
    std::cout << "Best value: " << bestValue << std::endl;
}

void TravellingSalesmanProblem::bruteForce_test() {
    if (array == NULL) {
        std::cout << "Array is empty!\n";
        return;
    }
    // przyjmijmy, �e zaczynamy od wierzcho�ka 0
    int source = 0;
    int bestValue = INT_MAX;
    int currentValue = 0;
    // najlepsza sciezka zapami�tana w tablicy wektorowej
    std::vector<int> bestTrip;
    // obecna sciezka zapamietana w t. wekt.
    std::vector<int> trip;

    // okreslenie domyslnej sciezki (wyjsciowej)
    for (int i = 0; i < V; i++) {
        trip.insert(trip.end(), i);
    }

    // okreslenie powrotu do punktu wyjscia dla domyslnej sciezki
    trip.insert(trip.end(), 0);

    // next_permutation z biblioteki <algorithm>, pozwala rearan�owa� elementy tab. wekt. tak aby uzyska� wszystkie mo�liwe jej permutacje
    while (std::next_permutation(trip.begin() + 1, trip.end() - 1)) {
        currentValue = 0;
        // dodawanie wag kraw�dzi grafu do ca�kowitego kosztu obecnej �cie�ki
        for (int i = 0; i < V; i++) {
            currentValue += array[trip[i]][trip[i + 1]];
        }

        // je�li obecna permutacja okaza�a si� korzystniejsza, to jest zapisywana jako owa
        if (currentValue < bestValue) {
            bestValue = currentValue;
            bestTrip = trip;
        }
    }
}

void TravellingSalesmanProblem::littleAlgorithm() {
    if (array == NULL) {
        std::cout << "Array is empty!\n";
        return;
    }
    // tablica wspolczynnikow standaryzacji
    std::vector<int> aFactor;
    std::vector<int> bFactor;

    int min = INT_MAX;

    int lowerBound = 0;

    for (int N = V; N > 0; N--) {
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                // szukamy minimum dla RZ�DU
                if (array[row][column] < min && row != column && array[row][column] != -1) {
                    min = array[row][column];
                }
            }
            // dodajemy minimum do tablicy wsp�cz. a
            aFactor.emplace_back(min);
            //std::cout << "a" << row << " " << min << std::endl;
            // ponowne ustalenie minimum
            min = INT_MAX;
        }

        // 1. krok do stworzenia C' - odj�cie od C wsp�czynnika a
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
            // dodajemy minimum do tablicy wsp�cz. a
            bFactor.emplace_back(min);
            //std::cout << "b" << column << " " << min << std::endl;
            // ponowne ustalenie minimum
            min = INT_MAX;
        }

        // 2. krok do stworzenia C' - odj�cie od C wsp�czynnika b
        for (int column = 0; column < V; column++) {
            for (int row = 0; row < V; row++) {
                // Cij - bi
                if (column != row && array[column][row] != -1) {
                    array[column][row] -= bFactor.at(row);
                }
            }
        }

        //print();

        // obliczenie dolnego oszacowania dla wszystkich rozwi�za�
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

        // stworzenie macierzy koszt�w rezygnacji dla tras "zerowych"
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
        //std::cout << "Resignation Matrix: " << std::endl;
        //for (int i = 0; i < V; i++) {
        //    for (int j = 0; j < V; j++) {
        //        std::cout << resignationArray[i][j] << " ";
        //    }
        //    std::cout << std::endl;
        //}

        // znalezienie d(kl) dla macierzy rezygnacji
        int d_kl = 0;
        // od razu zapisujemy row i column kt�re b�d� wykre�lane z macierzy
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

        // sprawdzenie warunku ko�cowego d_kl == 0
        if (d_kl == 0) {
            return;
        }

        // obliczenie dolnej granicy dla K2
        int lowerBound_K2 = lowerBound + d_kl;

        // tworzymy macierz zredukowan� C1
        // najpierw usuwamy rz�d
        for (int i = 0; i < V; i++) {
            array[k][i] = -1;
        }
        // potem kolumn�
        for (int j = 0; j < V; j++) {
            array[j][l] = -1;
        }
        // blokujemy podcykl tej samej �cie�ki
        array[l][k] = -1;

        //print();
    }

}

void TravellingSalesmanProblem::littleAlgorithm_test() {
    if (array == NULL) {
        std::cout << "Array is empty!\n";
        return;
    }
    // tablica wspolczynnikow standaryzacji
    std::vector<int> aFactor;
    std::vector<int> bFactor;

    int min = INT_MAX;

    int lowerBound = 0;

    for (int N = V; N > 0; N--) {
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                // szukamy minimum dla RZ�DU
                if (array[row][column] < min && row != column && array[row][column] != -1) {
                    min = array[row][column];
                }
            }
            // dodajemy minimum do tablicy wsp�cz. a
            aFactor.emplace_back(min);
            // ponowne ustalenie minimum
            min = INT_MAX;
        }

        // 1. krok do stworzenia C' - odj�cie od C wsp�czynnika a
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
            // dodajemy minimum do tablicy wsp�cz. a
            bFactor.emplace_back(min);
            // ponowne ustalenie minimum
            min = INT_MAX;
        }

        // 2. krok do stworzenia C' - odj�cie od C wsp�czynnika b
        for (int column = 0; column < V; column++) {
            for (int row = 0; row < V; row++) {
                // Cij - bi
                if (column != row && array[column][row] != -1) {
                    array[column][row] -= bFactor.at(row);
                }
            }
        }

        // obliczenie dolnego oszacowania dla wszystkich rozwi�za�
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

        // stworzenie macierzy koszt�w rezygnacji dla tras "zerowych"
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
        // od razu zapisujemy row i column kt�re b�d� wykre�lane z macierzy
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

        // sprawdzenie warunku ko�cowego d_kl == 0
        if (d_kl == 0) {
            return;
        }

        // obliczenie dolnej granicy dla K2
        int lowerBound_K2 = lowerBound + d_kl;

        // tworzymy macierz zredukowan� C1
        // najpierw usuwamy rz�d
        for (int i = 0; i < V; i++) {
            array[k][i] = -1;
        }
        // potem kolumn�
        for (int j = 0; j < V; j++) {
            array[j][l] = -1;
        }
        // blokujemy podcykl tej samej �cie�ki
        array[l][k] = -1;
    }
}


void TravellingSalesmanProblem::print()
{
    if (array == NULL) {
        std::cout << "Array is empty!\n";
        return;
    }

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

void TravellingSalesmanProblem::saveToFile()
{
    std::fstream file;
    file.open("generated.txt");

    file << V << "\n";

    for (int row = 0; row < V; row++) {
        for (int column = 0; column < V; column++) {
            file << array[row][column] << "    ";
        }
        file << "\n";
    }

    file.close();
}