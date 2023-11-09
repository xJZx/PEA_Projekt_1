#include "TravellingSalesmanProblem.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <queue>
#include <cmath>

TravellingSalesmanProblem::TravellingSalesmanProblem(int v, int low, int high)
{
    if (v == 0) {
        std::cout << "Matrix cannot be empty!\n";
        return;
    }

    this->V = v;
    // skierowany graf, w ktorym sciezki do obu tych samych wierzcholkow mog¹ byæ ró¿ne, w zale¿noœci z którego idziemy
    // czyli niesymetryczna macierz s¹siedztwa z wagami, dlatego taki wzór
    // ZMIENIÆ NA GRAF PE£NY!!!

    // wypelnienie macierzy maksymalnymi wartosciami
    matrix.resize(V, std::vector<int>(V, -1));

    int weight;
    for (int first = 0; first < V; first++)
    {
        for (int second = 0; second < V; second++)
        {
            if (first != second) {
                weight = (rand() % (high - low + 1)) + low;
                matrix[first][second] = weight;
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

        matrix.clear();
        matrix.resize(V, std::vector<int>(V));

        int weight;
        for (int first = 0; first < V; first++)
        {
            for (int second = 0; second < V; second++)
            {
                file >> weight;
                matrix[first][second] = weight;
                // bo graf asymetryczny skierowany
                //array[second][first] = weight;
            }      
        }
        file.close();
    }
    else
    {
        std::cout << "Error occurred!\n";
    }
}

TravellingSalesmanProblem::~TravellingSalesmanProblem() {
}

void TravellingSalesmanProblem::bruteForce() {
    if (matrix.empty() == true) {
        std::cout << "Array is empty!\n";
        return;
    }
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
            currentValue += matrix[trip[i]][trip[i + 1]];
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
    if (matrix.empty() == true) {
        std::cout << "Array is empty!\n";
        return;
    }
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
            currentValue += matrix[trip[i]][trip[i + 1]];
        }

        // jeœli obecna permutacja okaza³a siê korzystniejsza, to jest zapisywana jako owa
        if (currentValue < bestValue) {
            bestValue = currentValue;
            bestTrip = trip;
        }
    }
}

void TravellingSalesmanProblem::littleAlgorithm() {
    if (matrix.empty() == true) {
        std::cout << "Array is empty!\n";
        return;
    }
    
    // tablica odwiedzonych wierzcho³ków
    std::vector<int> visitedRow;
    std::vector<int> visitedColumn;

    // tablica wspolczynnikow standaryzacji
    std::vector<int> aFactor;
    std::vector<int> bFactor;

    // stworzenie kolejki
    std::queue<int> queueCost;
    std::queue<std::vector<std::vector<int>>> queueMatrix;

    // zaladowanie obecnej macierzy do kolejki
    queueMatrix.push(matrix);

    int min = INT_MAX;

    int lowerBound = INT_MAX;
    int lowerBound_K1 = 0;
    int iteration = 0;

    while (queueMatrix.empty() != true) {
        //for (int N = V; N > 0; N--) {
            for (int row = 0; row < V; row++) {
                for (int column = 0; column < V; column++) {
                    // szukamy minimum dla RZÊDU
                    if (matrix[row][column] < min && row != column && matrix[row][column] != -1) {
                        min = matrix[row][column];
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
                    if (row != column && matrix[row][column] != -1) {
                        matrix[row][column] -= aFactor.at(row);
                    }
                }
            }

            //print();

            for (int column = 0; column < V; column++) {
                for (int row = 0; row < V; row++) {
                    // szukamy minimum dla KOLUMNY
                    if (matrix[row][column] < min && row != column && matrix[row][column] != -1) {
                        min = matrix[row][column];
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
                    if (column != row && matrix[column][row] != -1) {
                        matrix[column][row] -= bFactor.at(row);
                    }
                }
            }

            print();

            // obliczenie dolnego oszacowania dla wszystkich rozwi¹zañ
            for (int i = 0; i < V; i++) {
                if (aFactor[i] == INT_MAX) {
                    lowerBound_K1 += 0;
                }
                else {
                    lowerBound_K1 += aFactor[i];
                }
                if (bFactor[i] == INT_MAX) {
                    lowerBound_K1 += 0;
                }
                else {
                    lowerBound_K1 += bFactor[i];
                }
            }

            //std::cout << "Lower bound on all tours: " << lowerBound << std::endl;

            // sprawdzenie czy wszystkie wspolczynniki sa rowne 0
            bool isFactorZero = false;
            for (int i = 0; i < aFactor.size(); i++) {
                if ((aFactor[i] == 0 || aFactor[i] == INT_MAX) && (bFactor[i] == 0 || bFactor[i] == INT_MAX)) {
                    isFactorZero = true;
                }
                else {
                    isFactorZero = false;
                    break;
                }
            }

            // sprawdzenie czy matrix posiada jedynie 0 i -1
            bool isMatrixZero = false;
            for (int i = 0; i < V; i++) {
                for (int j = 0; j < V; j++) {
                    if (matrix[i][j] == -1 || matrix[i][j] == 0) {
                        isMatrixZero = true;
                    }
                    else {
                        isMatrixZero = false;
                        break;
                    }
                } 
            }

            // wyczyszczenie tablic wektorowych
            aFactor.clear();
            bFactor.clear();

            // gdy wejdzie nowa macierz do kolejki
            checkNewMatrix:

            // stworzenie macierzy kosztów rezygnacji dla tras "zerowych"
            std::vector<std::vector<int>> resignationArray(V, std::vector<int>(V, -1));
            for (int row = 0; row < V; row++) {
                for (int column = 0; column < V; column++) {
                    if (matrix[row][column] == 0 && row != column && matrix[row][column] != -1) {
                        int minRow = INT_MAX;
                        int minColumn = INT_MAX;
                        for (int i = 0; i < V; i++) {
                            if (matrix[row][i] < minRow && i != column && matrix[row][i] != -1) {
                                minRow = matrix[row][i];
                            }
                        }
                        for (int j = 0; j < V; j++) {
                            if (matrix[j][column] < minColumn && j != row && matrix[j][column] != -1) {
                                minColumn = matrix[j][column];
                            }
                        }

                        // dodatkowe warunki
                        if ((minRow != INT_MAX && minColumn != INT_MAX) && (minRow != 0 || minColumn != 0)) {
                            resignationArray[row][column] = minRow + minColumn;
                        }
                        // jesli minRow == INT_MAX && minColumn != INT_MAX && minColumn!= 0
                        // kwestia czy "nieskonczonosc" + 0 powinno dawac 0
                        else if (minRow == INT_MAX && minColumn != INT_MAX /* && minColumn != 0*/) {
                            resignationArray[row][column] = minColumn;
                        }
                        // jesli minColumn == INT_MAX && minRow != INT_MAX && minRow!= 0
                        // kwestia czy "nieskonczonosc" + 0 powinno dawac 0
                        else if (minColumn == INT_MAX && minRow != INT_MAX /* && minRow != 0*/) {
                            resignationArray[row][column] = minRow;
                        }
                        // jesli zero w (k, l) ma w obu minimach zero, œwiadczy to ¿e jest ono podcyklem tamtych zer
                        else if (minRow == 0 && minColumn == 0) {
                            if (isFactorZero  && isMatrixZero) {
                                //matrix[row][column] = 0;
                                //matrix[column][row] = -1;
                                resignationArray[row][column] = -1;
                                //resignationArray[row][column] = 0;
                                //resignationArray[column][row] = -1;
                            }
                            else {
                                resignationArray[row][column] = 0;
                                //resignationArray[row][column] = -1;
                                //resignationArray[column][row] = 0;
                            }
                        }
                        /*else if (minRow == 0 && minColumn == 0) {
                            resignationArray[row][column] = -1;
                        }*/
                        // kwestia czy "nieskonczonosc" + 0 powinno dawac 0
                        /*else if ((minRow == INT_MAX && minColumn == 0) || (minRow == 0 && minColumn == INT_MAX)) {
                            resignationArray[row][column] = 0;
                        }*/
                        // jesli zero nie ma zadnego minimum, to jest ono ostatnim zerem, które zamyka nam cykl hamiltona
                        else if (minRow == INT_MAX && minColumn == INT_MAX) {
                            resignationArray[row][column] = 0;
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
            int d_kl = -1;
            // od razu zapisujemy row i column które bêd¹ wykreœlane z macierzy
            int k, l;
            for (int row = 0; row < V; row++) {
                for (int column = 0; column < V; column++) {
                    if (resignationArray[row][column] >= d_kl && resignationArray[row][column] != -1 && row != column) {
                        d_kl = resignationArray[row][column];
                        k = row;
                        l = column;
                    }
                }
            }
            std::cout << "d_kl = " << d_kl << std::endl;

            // sprawdzenie warunku koñcowego d_kl == -1
            if (d_kl == -1 || lowerBound_K1 > lowerBound) {
                // zapisanie nowego lowerBound
                if (lowerBound_K1 < lowerBound) {
                    lowerBound = lowerBound_K1;
                }
                // jesli siê oka¿e, ¿e lowerBound_K2 nie by³ mniejszy od lowerBound
                checkAgain:
                // usuniêcie z kolejki pierwszego elementu i sprawdzenie kolejnej macierzy
                queueMatrix.pop();
                if (queueMatrix.empty() == true) {
                    break;
                }
                // jesli kolejka kosztów jest mniejsza (lowerBound_K2 < lowerBound), to sprawdzamy kolejn¹ macierz   
                else if(queueCost.front() < lowerBound){
                    lowerBound_K1 = queueCost.front();
                    queueCost.pop();
                    matrix = queueMatrix.front();
                    visitedRow.push_back(-1);
                    visitedColumn.push_back(-1);
                    // i idziemy do nowej iteracji pêtli
                    if (iteration == 0) {
                        goto checkNewMatrix;
                    }
                    else {
                        continue;
                    }
                    //goto checkNewMatrix;
                }
                else if (queueCost.front() > lowerBound) {
                    queueCost.pop();
                    goto checkAgain;
                }
                else if (queueCost.front() == lowerBound) {
                    queueCost.pop();
                    matrix = queueMatrix.front();
                    visitedRow.push_back(-1);
                    visitedColumn.push_back(-1);
                    // i idziemy do nowej iteracji pêtli
                    if (iteration == 0) {
                        goto checkNewMatrix;
                    }
                    else {
                        continue;
                    }
                    //goto checkNewMatrix;
                }
            }
            else {
                visitedRow.push_back(k);
                visitedColumn.push_back(l);

                // obliczenie dolnej granicy dla K2
                int lowerBound_K2;
                if (iteration == 0) {
                    lowerBound_K2 = lowerBound_K1;
                }
                else {
                    lowerBound_K2 = lowerBound_K1 + d_kl;
                }

                std::vector<std::vector<int>> cloneArray(V, std::vector<int>(V, -1));

                // sprawdzenie czy mo¿e istnieæ optymalniejsza œcie¿ka
                // przepisanie obecnego arraya
                for (int row = 0; row < V; row++) {
                    for (int column = 0; column < V; column++) {
                        cloneArray[row][column] = matrix[row][column];
                    }
                }
                // usuwamy mozliwosc wybrania k, l dla K2
                cloneArray[k][l] = -1;

                // usuwamy podcykle dla obecnych tras
                int rightNeighbour, leftNeighbour;
                for (int v = 0; v < V; v++) {
                    for (int i = 0; i < visitedRow.size(); i++) {
                        if (visitedRow[i] == v) {
                            for (int j = 0; j < visitedColumn.size(); j++) {
                                if (visitedColumn[j] == v) {
                                    rightNeighbour = visitedColumn[i];
                                    leftNeighbour = visitedRow[j];

                                    cloneArray[rightNeighbour][leftNeighbour] = -1;
                                    cloneArray[leftNeighbour][rightNeighbour] = -1;
                                }
                            }
                        }
                    }
                }
                // wrzucamy do kolejki kolejn¹ macierz do sprawdzenia
                queueMatrix.push(cloneArray);
                // oraz odpowiadaj¹cy jej koszt
                queueCost.push(lowerBound_K2);
            }

            // tworzymy macierz zredukowan¹ C1
            // najpierw usuwamy rz¹d
            for (int i = 0; i < V; i++) {
                matrix[k][i] = -1;
            }
            // potem kolumnê
            for (int j = 0; j < V; j++) {
                matrix[j][l] = -1;
            }

            // usuwamy podcykle dla obecnych tras
            int rightNeighbour, leftNeighbour;
            for (int v = 0; v < V; v++) {
                for (int i = 0; i < visitedRow.size(); i++) {
                    if (visitedRow[i] == v) {
                        for (int j = 0; j < visitedColumn.size(); j++) {
                            if (visitedColumn[j] == v) {
                                rightNeighbour = visitedColumn[i];
                                leftNeighbour = visitedRow[j];

                                matrix[rightNeighbour][leftNeighbour] = -1;
                                matrix[leftNeighbour][rightNeighbour] = -1;
                            }
                        }
                    }
                }
            }

            // blokujemy podcykl tej samej œcie¿ki
            matrix[l][k] = -1;
           /* for (int i = 0; i < V; i++) {
                    matrix[i][k] = -1;
            }*/

            print();

            iteration++;
        //}
    }

    // algorytm na wytyczenie œcie¿ki z par szlaków miêdzy wierzcho³kami
    std::vector<int> path;
    int lastNeighbour = -1;
    for (int v = 0; v < V; v++) {
        for (int i = 0; i < visitedRow.size(); i++) {
            if (visitedRow[i] == -1 && visitedColumn[i] == -1) {
                break;
            }
            else {
                if (visitedRow[i] == 0 && v == 0) {
                    path.push_back(visitedRow[i]);
                    path.push_back(visitedColumn[i]);
                    lastNeighbour = visitedColumn[i];
                }
                else if (visitedRow[i] == lastNeighbour) {
                    path.push_back(visitedColumn[i]);
                    lastNeighbour = visitedColumn[i];
                    if (path[0] == lastNeighbour) {
                        break;
                    }
                }
            }   
        }
        if (path.empty() || lastNeighbour == 0) {
            break;
        }
    }

    std::cout << "The row path:    ";
    for (int i = 0; i < visitedRow.size(); i++) {
        std::cout << visitedRow[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "The column path: ";
    for (int i = 0; i < visitedColumn.size(); i++) {
        std::cout << visitedColumn[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "The path: ";
    for (int i = 0; i < path.size(); i++) {
        std::cout << path[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Total minimum cost: " << lowerBound << std::endl;
    return;

}

void TravellingSalesmanProblem::littleAlgorithm_test() {
    if (matrix.empty() == true) {
        std::cout << "Array is empty!\n";
        return;
    }

    // tablica odwiedzonych wierzcho³ków
    std::vector<int> visitedRow;
    std::vector<int> visitedColumn;

    // tablica wspolczynnikow standaryzacji
    std::vector<int> aFactor;
    std::vector<int> bFactor;

    int min = INT_MAX;

    int lowerBound = 0;

    for (int N = V; N > 0; N--) {
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                // szukamy minimum dla RZÊDU
                if (matrix[row][column] < min && row != column && matrix[row][column] != -1) {
                    min = matrix[row][column];
                }
            }
            // dodajemy minimum do tablicy wspó³cz. a
            aFactor.emplace_back(min);
            //std::cout << "a" << row << " " << min << std::endl;
            // ponowne ustalenie minimum
            min = INT_MAX;
        }

        // 1. krok do stworzenia C' - odjêcie od C wspó³czynnika a
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                // Cij - ai
                if (row != column && matrix[row][column] != -1) {
                    matrix[row][column] -= aFactor.at(row);
                }
            }
        }

        for (int column = 0; column < V; column++) {
            for (int row = 0; row < V; row++) {
                // szukamy minimum dla KOLUMNY
                if (matrix[row][column] < min && row != column && matrix[row][column] != -1) {
                    min = matrix[row][column];
                }
            }
            // dodajemy minimum do tablicy wspó³cz. a
            bFactor.emplace_back(min);
            //std::cout << "b" << column << " " << min << std::endl;
            // ponowne ustalenie minimum
            min = INT_MAX;
        }

        // 2. krok do stworzenia C' - odjêcie od C wspó³czynnika b
        for (int column = 0; column < V; column++) {
            for (int row = 0; row < V; row++) {
                // Cij - bi
                if (column != row && matrix[column][row] != -1) {
                    matrix[column][row] -= bFactor.at(row);
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
                if (matrix[row][column] == 0 && row != column && matrix[row][column] != -1) {
                    int minRow = INT_MAX;
                    int minColumn = INT_MAX;
                    for (int i = 0; i < V; i++) {
                        if (matrix[row][i] < minRow && i != column && matrix[row][i] != -1) {
                            minRow = matrix[row][i];
                        }
                    }
                    for (int j = 0; j < V; j++) {
                        if (matrix[j][column] < minColumn && j != row && matrix[j][column] != -1) {
                            minColumn = matrix[j][column];
                        }
                    }
                    if ((minRow != INT_MAX && minColumn != INT_MAX) && (minRow != 0 || minColumn != 0)) {
                        resignationArray[row][column] = minRow + minColumn;
                    }
                    // jesli minRow == INT_MAX && minColumn != INT_MAX && minColumn!= 0
                    else if (minRow == INT_MAX && minColumn != INT_MAX && minColumn != 0) {
                        resignationArray[row][column] = minColumn;
                    }
                    // jesli minColumn == INT_MAX && minRow != INT_MAX && minRow!= 0
                    else if (minColumn == INT_MAX && minRow != INT_MAX && minRow != 0) {
                        resignationArray[row][column] = minRow;
                    }
                    // jesli zero w (k, l) ma w obu minimach zero, œwiadczy to ¿e jest ono podcyklem tamtych zer
                    else if (minRow == 0 && minColumn == 0) {
                        resignationArray[row][column] = -1;
                    }
                    // kwestia czy "nieskonczonosc" + 0 powinno dawac 0
                    else if ((minRow == INT_MAX && minColumn == 0) || (minRow == 0 && minColumn == INT_MAX)) {
                        resignationArray[row][column] = 0;
                    }
                    // jesli zero nie ma zadnego minimum, to jest ono ostatnim zerem, które zamyka nam cykl hamiltona
                    else if (minRow == INT_MAX && minColumn == INT_MAX) {
                        resignationArray[row][column] = 0;
                    }
                }
            }
        }

        // znalezienie d(kl) dla macierzy rezygnacji
        int d_kl = -1;
        // od razu zapisujemy row i column które bêd¹ wykreœlane z macierzy
        int k, l;
        for (int row = 0; row < V; row++) {
            for (int column = 0; column < V; column++) {
                if (resignationArray[row][column] >= d_kl && resignationArray[row][column] != -1 && row != column) {
                    d_kl = resignationArray[row][column];
                    k = row;
                    l = column;
                }
            }
        }

        visitedRow.push_back(k);
        visitedColumn.push_back(l);

        // sprawdzenie warunku koñcowego d_kl == 0
        if (d_kl == -1) {
            return;
        }

        // obliczenie dolnej granicy dla K2
        int lowerBound_K2 = lowerBound + d_kl;

        // tworzymy macierz zredukowan¹ C1
        // najpierw usuwamy rz¹d
        for (int i = 0; i < V; i++) {
            matrix[k][i] = -1;
        }
        // potem kolumnê
        for (int j = 0; j < V; j++) {
            matrix[j][l] = -1;
        }

        // usuwamy podcykle dla obecnych tras
        int rightNeighbour, leftNeighbour;
        for (int v = 0; v < V; v++) {
            for (int i = 0; i < visitedRow.size(); i++) {
                if (visitedRow[i] == v) {
                    for (int j = 0; j < visitedColumn.size(); j++) {
                        if (visitedColumn[j] == v) {
                            rightNeighbour = visitedColumn[i];
                            leftNeighbour = visitedRow[j];

                            matrix[rightNeighbour][leftNeighbour] = -1;
                            matrix[leftNeighbour][rightNeighbour] = -1;
                        }
                    }
                }
            }
        }

        // blokujemy podcykl tej samej œcie¿ki
        matrix[l][k] = -1;
    }
}

// mask, to iloœæ bitów w jakiej znajduje siê miasto
// pos, to obecne miasto w jakim siê znajdujemy
void TravellingSalesmanProblem::dynamicProgramming()
{
    // tablica 2^n, V, wype³nienie jej -1
    dp.clear();
    dp.resize(1 << V, std::vector<int>(V, -1));
    // tablica dla œcie¿ki
    parents.clear();
    parents.resize(1 << V, std::vector<int>(V, -1));
    //parents.push_back(0);

    std::cout << "Cost of the path: " << dynamicProgrammingRecursion(1, 0) << std::endl;

    // startujemy odwzorowywaæ œcie¿kê w tym samym punkcie startowym, w którym zaczêliœmy algorytm
    int mask = 1;
    int pos = 0;
    std::vector<int> path;
    path.push_back(0);

    while (mask) {
        // jeœli obecny rodzic bêdzie równy -1, to znaczy ¿e wyznaczyliœmy ju¿ ca³¹ œcie¿kê
        if (parents[mask][pos] != -1){
            int prev_pos = parents[mask][pos];
            path.push_back(prev_pos);
            mask = mask | (1 << prev_pos);
            pos = prev_pos;
        }
        else {
            break;
        } 
    }

    path.push_back(0);

    std::cout << "Best trip: ";
    for (int i = 0; i < path.size(); i++) {
        std::cout << path[i] << " ";
    }
    std::cout << std::endl;
    
}

int TravellingSalesmanProblem::dynamicProgrammingRecursion(int mask, int pos)
{
    // sprawdzenie czy wszystkie miasta by³y ju¿ odwiedzone
    // 1 << V, to inaczej 2^V (operacje bitowe)
    if (mask == (1 << V) - 1) {
        return matrix[pos][0];
    }

    // sprawdzenie czy problem zostal ju¿ rozwi¹zany
    if (dp[mask][pos] != -1) {
        return dp[mask][pos];
    }

    int minCost = INT_MAX;

    // odwiedzenie nieodwiedzonych miast i znalezienie najkrótszej œcie¿ki
    for (int city = 0; city < V; city++) {
        if ((mask & (1 << city)) == 0) {
            int newMinCost = matrix[pos][city] + dynamicProgrammingRecursion(mask | (1 << city), city);
            //minCost = std::min(minCost, newMinCost);
            if (newMinCost < minCost) {
                minCost = newMinCost;
                //pathDynamic.push_back(city);
                parents[mask][pos] = city;
            }
        }
    }

    dp[mask][pos] = minCost;

    return minCost;
}

void TravellingSalesmanProblem::dynamicProgramming_test() {
    // tablica 2^n, V, wype³nienie jej -1
    dp.clear();
    dp.resize(1 << V, std::vector<int>(V, -1));
    // tablica dla œcie¿ki
    parents.clear();
    parents.resize(1 << V, std::vector<int>(V, -1));
    //parents.push_back(0);

    dynamicProgrammingRecursion(1, 0);

    // startujemy odwzorowywaæ œcie¿kê w tym samym punkcie startowym, w którym zaczêliœmy algorytm
    int mask = 1;
    int pos = 0;
    std::vector<int> path;
    path.push_back(0);

    while (mask) {
        // jeœli obecny rodzic bêdzie równy -1, to znaczy ¿e wyznaczyliœmy ju¿ ca³¹ œcie¿kê
        if (parents[mask][pos] != -1) {
            int prev_pos = parents[mask][pos];
            path.push_back(prev_pos);
            mask = mask | (1 << prev_pos);
            pos = prev_pos;
        }
        else {
            break;
        }
    }

    path.push_back(0);
}


void TravellingSalesmanProblem::print()
{
    if (matrix.empty() == true) {
        std::cout << "Array is empty!\n";
        return;
    }

    std::cout << "Graph stored in the matrix: \n";
    for (int first = 0; first < V; first++)
    {
        for (int second = 0; second < V; second++)
        {
            std::cout << std::setw(4) << matrix[first][second];
        }

        std::cout << "\n";
    }
    std::cout << std::endl;
}

void TravellingSalesmanProblem::saveToFile()
{
    std::ofstream file;
    file.open("generated.txt");

    file << V << "\n";

    for (int row = 0; row < V; row++) {
        for (int column = 0; column < V; column++) {
            file << matrix[row][column] << "    ";
        }
        file << "\n";
    }

    file.close();
}