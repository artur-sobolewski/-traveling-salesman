#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;

int minimalna = INT_MAX;
int optimum;
int comb = 1;
int UPPERBOUND;
int ile = 0;
int close = 0;
int ilosc = 0;

int* nearest_neigh(int** edge_matrix, int num_of_cities)
{
    int seq = 2;
    int length;
    int act_node = 0;
    int min_edge = INT_MAX;
    int neighbour = 0;
    int* seqence = new int[num_of_cities];
    int* check = new int[num_of_cities];
    for (int i = 0; i < num_of_cities; i++)
    {
        check[i] = 0;
    }
    length = 0;
    while (seq <= num_of_cities)
    {
        act_node = neighbour;
        min_edge = INT_MAX;
        check[0] = 1;
        for (int j = 0; j < num_of_cities; j++)
        {
            // Czy długość krawędzi mniejsza niż aktualnie najmniejsza
            // Czy dany wierzchołek nie został rozpatrzony
            // Czy dany wierzchołek nie jest wierzchołkiem, z którego krawędź wychodzi
            // Czy istnieje krawędź pomiędzy wierzchołkami
            if ((edge_matrix[act_node][j] < min_edge) && (check[j] == 0) && (act_node != j) && (edge_matrix[act_node][j] != 0))
            {
                min_edge = edge_matrix[act_node][j];
                neighbour = j;
            }
        }
        check[neighbour] = seq;
        seq++;
    }
    for (int i = 1; i < num_of_cities + 1; i++)
    {
        for (int j = 0; j < num_of_cities; j++)
        {
            if (check[j] == i)
                seqence[i - 1] = j;
        }
    }
    return seqence;
}

void copySeq(int* seq, int* directSeq, int size)
{
    for (int i = 0; i < size; i++)
    {
        directSeq[i] = seq[i];
    }
}

void insert(int* sequence, int i, int j, int size)
{

    int pom = sequence[j];
    for (int it = j; it > i; it--)
    {
        sequence[it] = sequence[it - 1];
    }
    sequence[i] = pom;
}
void swap(int* sequence, int i, int j, int size)
{
    int pom = sequence[i];
    sequence[i] = sequence[j];
    sequence[j] = pom;
}
void invert(int* sequence, int i, int j, int size)
{
    int cor = 0;
    if ((j - i) % 2 != 0)
    {
        cor = 1;
    }
    for (int it = 0; it < ((j - i) / 2) + cor; it++)
    {
        swap(sequence, i + it, j - it, size);
    }
}



void displayPath(int* seq, int size)
{
    for (int i = 0; i < size - 1; i++)
    {
        cout << seq[i] << " -> ";
    }
    cout << seq[size - 1] << endl;
}

struct Individuals
{
    int* chromosome;
    float value;
    int length;
    float percent1;
    float percent2;
};


void sortowanie_babelkowe(Individuals tab[], int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 1; j < n - i; j++)
            if (tab[j - 1].length > tab[j].length)
                swap(tab[j - 1], tab[j]);
}

class GenethicAlgorith
{
public:
    int k1, k2, pom, pomx, index;
    float pc, pm;
    float prob;
    int size, pomocnicza;
    int populationSize;
    int** edge_matrix;
    Individuals* population;
    Individuals best;
    Individuals* newAndOldPop;
    int* mappingTab;
    bool* ifChecked1;
    bool* ifChecked2;
    int* nearestNei;
    int nearestLen;
    Individuals child1;
    Individuals child2;
    GenethicAlgorith(int size, int** matrix, int popSize, float pc, float pm);
    void generatePerm(int* perm);
    float setValue(Individuals& ind);
    void startGA(int iteracje);
    void displayPath(int* perm);
    void PMX(int* ind1, int* ind2, int* child1, int* child2);
    int calculateLength(int* perm);
};

GenethicAlgorith::GenethicAlgorith(int size, int** matrix, int popSize, float pc, float pm)
{
    this->size = size;
    edge_matrix = matrix;
    populationSize = popSize;
    this->pc = pc;
    this->pm = pm;
    newAndOldPop = new Individuals[popSize * 2];
    mappingTab = new int[size + 1];
    best.chromosome = new int[size];
    ifChecked1 = new bool[size];
    ifChecked2 = new bool[size];
    for (int i = 0; i < size; i++)
    {
        ifChecked1[i] = false;
        ifChecked2[i] = false;
    }
    for (int i = 0; i < populationSize * 2; i++)
    {
        newAndOldPop[i].chromosome = new int[size];
        for (int j = 0; j < size; j++)
            newAndOldPop[i].chromosome[j] = j;
    }
    nearestNei = nearest_neigh(edge_matrix, size);
    nearestLen = calculateLength(nearestNei);
}

int GenethicAlgorith::calculateLength(int* perm)
{
    int len = 0;
    for (int i = 0; i < size; i++)
    {
        len += edge_matrix[perm[i]][perm[(i + 1) % size]];
    }
    return len;
}

void GenethicAlgorith::generatePerm(int* perm)
{
    for (int i = 0; i < size; i++)
    {
        perm[i] = i;
    }
    int j;
    for (int i = 0; i < size; i++)
    {
        j = rand() % (size);
        swap(perm, i, j, size);
    }
}

float GenethicAlgorith::setValue(Individuals& ind)
{
    int len = 0;

    for (int i = 0; i < size - 1; i++)
    {
        len += edge_matrix[ind.chromosome[i]][ind.chromosome[i + 1]];
    }
    len += edge_matrix[ind.chromosome[size - 1]][ind.chromosome[0]];
    ind.length = len;
    return (float)nearestLen / len;
}

void GenethicAlgorith::PMX(int* ind1, int* ind2, int* child1, int* child2)
{
    int x;
    do
    {
        k1 = rand() % size;
        k2 = rand() % size;
    } while (k1 == k2);
    if (k1 > k2)
    {
        pom = k1;
        k1 = k2;
        k2 = pom;
    }
    for (int i = 0; i <= k2 - k1; i++)
    {
        mappingTab[ind1[k1 + i]] = ind2[k1 + i];

        ifChecked2[ind1[k1 + i]] = true;
        ifChecked1[ind2[k1 + i]] = true;

        child1[k1 + i] = ind2[k1 + i];
        child2[k1 + i] = ind1[k1 + i];
    }
    for (int i = 0; i < size; i++)
    {
        if (i < k1 || i > k2)
        {
            if (ifChecked2[ind2[i]] == true)
            {
                x = mappingTab[ind2[i]];
                while (ifChecked2[x] == true)
                {
                    x = mappingTab[x];
                }
                child2[i] = x;
            }
            else
            {
                child2[i] = ind2[i];
            }
            if (ifChecked1[ind1[i]] == true)
            {
                x = ind1[i];
                do
                {
                    index = 0;
                    while (mappingTab[index] != x && index < size)
                    {
                        index++;
                    }
                    x = index;
                } while (ifChecked1[x] == true);
                child1[i] = x;
            }
            else
            {
                child1[i] = ind1[i];
            }
        }
    }
    for (int i = 0; i < size; i++)
    {
        ifChecked1[i] = false;
        ifChecked2[i] = false;
        mappingTab[i] = -1;
    }
    mappingTab[size] = -1;
    do
    {
        k1 = rand() % size;
        k2 = rand() % size;
    } while (k1 == k2);
    if (k1 > k2)
    {
        pom = k1;
        k1 = k2;
        k2 = pom;
    }
    
    pomocnicza = rand() % 1000;
    prob = (float)pomocnicza / 1000;
    if (prob < pm)
    {
        insert(child1, k1, k2, size);
    }
    pomocnicza = rand() % 1000;
    prob = (float)pomocnicza / 1000;
    if (prob < pm)
    {
        insert(child2, k1, k2, size);
    }
}


void GenethicAlgorith::startGA(int iteracje)
{
    bool ok;
    int pos = 0;
    int minLen = INT_MAX;
    Individuals crossInd1, crossInd2;
    int x;
    int k_prev = -1, k = 0;
    int ile = 0;
    float prevPercent = 0, percent, total = 0;
    int crossNum = 0;
    for (int i = 0; i < populationSize; i++)
    {
        generatePerm(newAndOldPop[i].chromosome);
        newAndOldPop[i].value = setValue(newAndOldPop[i]);
        if (newAndOldPop[i].length < minLen)
        {
            minLen = newAndOldPop[i].length;
            best = newAndOldPop[i];
            cout << "START " << minLen << " (" << (float)minLen / optimum * 100 << "%)" << endl;
        }
    }
    
    for (int i = 0; i < iteracje; i++)
    {
        //Selekcja
        total = 0;
        for (int j = 0; j < populationSize; j++)
        {
            newAndOldPop[j].value = setValue(newAndOldPop[j]);
            total += newAndOldPop[j].value;
        }

        prevPercent = 0;
        for (int j = 0; j < populationSize; j++)
        {
            newAndOldPop[j].percent1 = prevPercent;
            percent = newAndOldPop[j].value / total * 1000.0;
            newAndOldPop[j].percent2 = prevPercent + percent;
            prevPercent += percent;
        }
        ile = 0;
        crossNum = 0;
        ok = false;
        while (crossNum < populationSize)
        {
            x = rand() % 10000;
            percent = (float)x / 100;

            for (int k = 0; k < populationSize; k++)
            {
                if ((percent >= newAndOldPop[k].percent1) && (percent < newAndOldPop[k].percent2))
                {
                    if (k == k_prev)
                    {
                        break;
                    }
                    else
                    {
                        if (ile == 0)
                        {
                            k_prev = k;
                            crossInd1 = newAndOldPop[k];
                            ile++;
                        }
                        else if (ile == 1)
                        {
                            crossInd2 = newAndOldPop[k];
                            PMX(crossInd1.chromosome, crossInd2.chromosome,
                                newAndOldPop[populationSize + crossNum].chromosome, newAndOldPop[populationSize + crossNum + 1].chromosome);
                            k_prev = -1;
                            crossNum += 2;
                            ile = 0;
                            ok = false;
                        }
                    }
                }
            }
        }

        //Po krzyżowaniu
        total = 0;
        for (int j = 0; j < populationSize * 2; j++)
        {
            newAndOldPop[j].value = setValue(newAndOldPop[j]);
            total += newAndOldPop[j].value;
            if (newAndOldPop[j].length < minLen)
            {
                minLen = newAndOldPop[j].length;
                cout << i << " Dlugosc: " << minLen << " (" << (float)minLen / optimum * 100 << "%)" << endl;
                best = newAndOldPop[j];
            }
        }

        sortowanie_babelkowe(newAndOldPop, populationSize * 2);
       
    }
}



void GenethicAlgorith::displayPath(int* perm)
{
    cout << endl;
    for (int i = 0; i <= size; i++)
    {
        cout << perm[i % size] << " ";
    }
    cout << endl;
}


enum Operacja
{
    DEFAULT,
    INSERT,
    INVERT,
    SWAP
};
struct Atrybut
{
    Operacja operation;
    int i;
    int j;
    int kadencja;
    void toSting()
    {
        cout << operation << " " << i << " " << j << " " << kadencja << endl;
    }
};

class TabuSearch
{
public:
    int size;
    int** edge_matrix;
    int* new_path;
    int* act_path;
    int* final_path;
    int lastPos;
    int* pom_seq_insert;
    int* pom_seq_swap;
    int* pom_seq_invert;
    int* perm;
    bool swapDozwolone;
    bool insertDozwolone;
    bool invertDozwolone;
    Atrybut* Tabu;
    TabuSearch(int size, int** matrix);
    void removeFromTabu(int index, int last);
    int calculateCost(int* seq);
    void startTS(int itaration);
    int* generatePerm();
    int* bestNeig(int* seq, int act_len, int* new_len, int i, int j, int* lastOnTabu);
};

TabuSearch::TabuSearch(int size, int** matrix)
{
    this->size = size;
    edge_matrix = matrix;
    new_path = new int[size];
    act_path = new int[size];
    final_path = new int[size];
    Tabu = new Atrybut[(size) * (size - 1) * 3];
    lastPos = 0;
    pom_seq_insert = new int[size];
    pom_seq_invert = new int[size];
    pom_seq_swap = new int[size];
    for (int i = 0; i < (size) * (size - 1) * 3; i++)
    {
        Tabu[i] = { Operacja::DEFAULT, -1, -1, -1 };
    }
    perm = new int[size];
    for (int i = 0; i < size; i++)
    {
        perm[i] = i;
    }
}

int* TabuSearch::generatePerm()
{
    for (int i = 0; i <= size - 2; i++)
    {
        perm[i] = i;
    }
    int j;
    for (int i = 0; i <= size - 2; i++)
    {
        j = rand() % (size);
        swap(perm, i, j, size);
    }
    return perm;
}

int TabuSearch::calculateCost(int* seq)
{
    int len = 0;
    for (int i = 0; i < size; i++)
    {
        len += edge_matrix[seq[i]][seq[(i + 1) % size]];
    }
    return len;
}

int* TabuSearch::bestNeig(int* seq, int act_len, int* new_len, int i_ind, int j_ind, int* lastOnTabu)
{
    int wybor = 0;
    copySeq(seq, pom_seq_swap, size);
    swapDozwolone = true;
    int it = 0;
    int len;
    int min = INT_MAX;
    Atrybut atrybut = { Operacja::DEFAULT, i_ind, j_ind, (size / 2) + (rand() % (size / 2)) };

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (i != j)
            {
                copySeq(seq, pom_seq_swap, size);
                swap(pom_seq_swap, i, j, size);
                len = calculateCost(pom_seq_swap);
                if (len < min)
                {
                    min = len;
                    atrybut = { Operacja::SWAP, i, j, (size / 2) + (rand() % (size / 2)) };
                }
            }
        }
    }

    while (it < *lastOnTabu)
    {
        if (Tabu[it].i == atrybut.i && Tabu[it].j == atrybut.j)
        {
            if (Tabu[it].operation == Operacja::SWAP)
            {
                swapDozwolone = false;
            }
        }
        if (swapDozwolone == true)
        {
            it++;
        }
        else
        {
            break;
        }
    }
    copySeq(seq, pom_seq_swap, size);
    swap(pom_seq_swap, atrybut.i, atrybut.j, size);
    len = calculateCost(pom_seq_swap);
    if (swapDozwolone == false)
    {
        if (len < act_len)
        {
            *new_len = len;
            return pom_seq_swap;
        }
    }
    else
    {
        *new_len = len;
        if (swapDozwolone == true)
        {
            Tabu[*lastOnTabu] = atrybut;
            *lastOnTabu = *lastOnTabu + 1;
        }
        return pom_seq_swap;
    }
}

void TabuSearch::removeFromTabu(int index, int last)
{
    int i = index;
    while (i < last)
    {
        Tabu[i] = Tabu[i + 1];
        i++;

        if (i == (size) * (size - 1) * 3 - 1)
        {
            Tabu[i] = { Operacja::DEFAULT, -1, -1, -1 };
            return;
        }
    }
}



void TabuSearch::startTS(int it_param)
{
    int MAX_iteration = size * it_param;
    int zam;
    int i, j;
    int freeOnTabu = 0;
    int* wsk_freeOnTabu = &freeOnTabu;
    int iteration = 0;
    act_path = nearest_neigh(edge_matrix, size);
    int act_len = calculateCost(act_path);
    int new_len = 0;
    int min_len = act_len;
    int* wsk_new_len = &new_len;
    copySeq(act_path, final_path, size);
    cout << "START  " << min_len << " (" << (float)min_len / optimum * 100 << "%)" << endl;
    int correctTabuIterator = 0;
    while (iteration < MAX_iteration)
    {

        //new_path = bestNeig(act_path, min_len, wsk_new_len, wsk_freeOnTabu);
//           do
//           {
        do
        {
            i = rand() % size;
            j = rand() % size;
            if (i > j)
            {
                zam = i;
                i = j;
                j = zam;
            }
        } while (i == j);
        new_path = bestNeig(act_path, min_len, wsk_new_len, i, j, wsk_freeOnTabu);
        //            } while (act_len == new_len);
        if (new_len < act_len)
        {
            copySeq(new_path, act_path, size);
            act_len = new_len;
        }


        if (new_len < min_len)
        {
            min_len = new_len;
            copySeq(act_path, final_path, size);
            cout << iteration << "  " << min_len << " (" << (float)min_len / optimum * 100 << "%)" << endl;
            ile = 0;
        }
        else
        {
            ile++;
        }

        while (correctTabuIterator < freeOnTabu)
        {
            Tabu[correctTabuIterator].kadencja--;
            if (Tabu[correctTabuIterator].kadencja < 1)
            {
                removeFromTabu(correctTabuIterator, freeOnTabu);
                if (freeOnTabu != 0)
                {
                    freeOnTabu--;
                }
            }
            else
            {
                correctTabuIterator++;
            }
        }
        correctTabuIterator = 0;
        iteration++;
    }
}


class SimulatedAnnealing
{
public:
    int size;
    int** edge_matrix;
    int* sequence;
    int* act_path;
    int* final_path;
    int* pom_seq_insert;
    int* pom_seq_swap;
    int* pom_seq_invert;
    int* perm;
    SimulatedAnnealing(int size, int** matrix);
    int* bestNeig(int* seq, int* len, int i, int j);
    int calculateCost(int* seq);
    void startSA(double start_temp, double cooling_param, double min_temp);
    int* generatePerm();
};

SimulatedAnnealing::SimulatedAnnealing(int size, int** matrix)
{
    this->size = size;
    edge_matrix = matrix;
    sequence = new int[size];
    act_path = new int[size];
    final_path = new int[size];
    pom_seq_insert = new int[size];
    pom_seq_invert = new int[size];
    pom_seq_swap = new int[size];
    perm = new int[size];
    for (int i = 0; i < size; i++)
    {
        perm[i] = i;
    }
}



int* SimulatedAnnealing::generatePerm()
{
    for (int i = 0; i <= size - 2; i++)
    {
        perm[i] = i;
    }
    int j;
    for (int i = 0; i <= size - 2; i++)
    {
        j = rand() % (size - 1);
        swap(perm, i, j, size);
    }
    return perm;
}

int* SimulatedAnnealing::bestNeig(int* seq, int* len, int i, int j)
{
    int wybor = 0;
    copySeq(seq, pom_seq_insert, size);
    copySeq(seq, pom_seq_swap, size);
    copySeq(seq, pom_seq_invert, size);
    insert(pom_seq_insert, i, j, size);
    int insertLen = calculateCost(pom_seq_insert);
    int min = insertLen;
    swap(pom_seq_swap, i, j, size);
    int swapLen = calculateCost(pom_seq_swap);
    if (swapLen < insertLen)
    {
        min = swapLen;
        wybor = 1;
    }
    invert(pom_seq_invert, i, j, size);
    int invertLen = calculateCost(pom_seq_invert);
    if (invertLen < min)
    {
        min = invertLen;
        wybor = 2;
    }
    if (wybor == 0)
    {
        *len = insertLen;
        return pom_seq_insert;
    }
    else if (wybor == 1)
    {
        *len = swapLen;
        return pom_seq_swap;
    }
    else if (wybor == 2)
    {
        *len = invertLen;
        return pom_seq_invert;
    }
}

int SimulatedAnnealing::calculateCost(int* seq)
{
    int len = 0;
    for (int i = 0; i < size; i++)
    {
        len += edge_matrix[seq[i]][seq[(i + 1) % size]];
    }
    return len;
}

void SimulatedAnnealing::startSA(double start_temp, double cooling_param, double min_temp)
{
    int zam;
    double temp = start_temp;
    act_path = nearest_neigh(edge_matrix, size);
    copySeq(act_path, final_path, size);
    int new_len = calculateCost(act_path);
    int* wsk_new_len = &new_len;
    int min_len = new_len;
    int prev_len = new_len;
    int i, j, iteracja = 0, delta_len, ile = 0;
    double p, r;
    cout << "START  " << new_len << " (" << (float)new_len / optimum * 100 << "%)" << endl;
    while (temp > min_temp)
    {

        for (int inner_it = 0; inner_it < size; inner_it++)
        {
            iteracja++;

            do
            {
                i = rand() % size;
                j = rand() % size;
                if (i > j)
                {
                    zam = i;
                    i = j;
                    j = zam;
                }

            } while (i == j);

            if (ile > 1000)
            {
                sequence = generatePerm();
                new_len = calculateCost(sequence);
                ile = 0;
            }
            else
            {
                sequence = bestNeig(act_path, wsk_new_len, i, j);
            }
            delta_len = new_len - prev_len;
            if (delta_len < 0)
            {
                prev_len = new_len;
                copySeq(sequence, act_path, size);
            }
            else
            {
                ile++;
                p = exp(-(delta_len / temp));
                r = (float)(rand() % 1000) / 1000;
                if (r < p)
                {
                    prev_len = new_len;
                    copySeq(sequence, act_path, size);
                }
            }
            if (new_len < min_len)
            {
                min_len = new_len;
                copySeq(act_path, final_path, size);
                cout << iteracja << "  " << min_len << " (" << (float)min_len / optimum * 100 << "%)" << endl;
                ile = 0;
            }
        }
        temp = cooling_param * temp;
    }
    return;
}




void bruteForce(int** edge_matrix, int* permutacje, int N, int size, int* seq)
{
    int min;
    int pom;
    if (N == 0)
    {
        float prd = 0;
        int len = 0;
        for (int i = 0; i < size; i++)
        {
            len += edge_matrix[permutacje[i]][permutacje[(i + 1) % size]];
        }
        if (len < minimalna)
        {
            minimalna = len;
            prd = 100.0 * (float)len / (float)optimum;
            cout << len << " " << prd << "%" << endl;
            for (int j = 0; j < size; j++)
            {
                seq[j] = permutacje[j];
            }
        }
    }
    else
        for (int i = 0; i < N; i++)
        {
            bruteForce(edge_matrix, permutacje, N - 1, size, seq);
            if (N % 2 == 1)
                swap(permutacje[0], permutacje[N - 1]);
            else
                swap(permutacje[i], permutacje[N - 1]);
        }
    return;
}

void startBruteF(int** edge_matrix, int* num_of_cities, int* seq)
{
    int* permutacje = new int[*num_of_cities];
    for (int i = 0; i < *num_of_cities; i++)
    {
        permutacje[i] = i;
    }
    bruteForce(edge_matrix, permutacje, *num_of_cities - 1, *num_of_cities, seq);

    return;
}


class BandBClass
{
public:
    int* sequence;
    int** edge_matrix;
    bool* checked;
    int** kolejka;
    int* path;
    int size;
    int minimum = INT_MAX;
    BandBClass(int size, int** matrix);
    void displayPath(void);
    void b_and_b_start();
    void b_and_b_loop(int N, int length, int lowerBound);
};
BandBClass::BandBClass(int size, int** matrix)
{
    this->size = size;
    this->checked = new bool[size];
    this->sequence = new int[size];
    this->edge_matrix = matrix;
    this->kolejka = new int* [size];
    this->path = new int[size];
    for (int i = 0; i < size; i++)
        sequence[i] = 0;
    for (int i = 0; i < size; i++)
        checked[i] = false;
    for (int i = 0; i < size; i++)
    {
        kolejka[i] = new int[2];
    }
}

void BandBClass::displayPath()
{
    for (int i = 0; i < size; i++)
        cout << path[i] << " -> ";
    cout << path[0] << endl;
}

void BandBClass::b_and_b_loop(int N, int length, int lowerBound)
{
    int len;
    int actLowerBound;
    if (N == size)
    {
        len = length;
        len += edge_matrix[sequence[N - 1]][sequence[0]];
        if (len < minimum)
        {
            minimum = len;
            cout << len << " " << 100.0 * len / optimum << "%" << endl;
            for (int i = 0; i < size; i++)
            {
                path[i] = sequence[i];
            }
        }
        return;
    }
    else
    {
        for (int j = 0; j < size; j++)
        {
            if (checked[j] != true)
            {
                len = length;
                actLowerBound = lowerBound;
                len += edge_matrix[sequence[N - 1]][j];
                if (N == 1)
                {
                    actLowerBound -= (kolejka[N - 1][0] + kolejka[j][1]) / 2;
                }
                else
                {
                    actLowerBound -= (kolejka[N - 1][1] + kolejka[j][1]) / 2;
                }
                if (actLowerBound + len < minimum)
                {
                    sequence[N] = j;
                    checked[j] = true;
                    b_and_b_loop(N + 1, len, actLowerBound);

                }
                checked[j] = false;
            }
        }
    }
}

void BandBClass::b_and_b_start()
{
    int lowerBound = 0;
    for (int i = 0; i < size; i++)
    {
        int firstMin = INT_MAX, secondMin = INT_MAX;
        for (int j = 0; j < size; j++)
        {
            if (edge_matrix[i][j] < firstMin && i != j)
            {
                secondMin = firstMin;
                firstMin = edge_matrix[i][j];
            }
            else if (edge_matrix[i][j] < secondMin && i != j)
            {
                secondMin = edge_matrix[i][j];
            }
        }
        kolejka[i][0] = firstMin;
        kolejka[i][1] = secondMin;
        lowerBound += firstMin + secondMin;
    }
    if (lowerBound % 2 == 1)
        lowerBound = lowerBound / 2 + 1;
    else lowerBound = lowerBound / 2;
    checked[0] = true;
    b_and_b_loop(1, 0, lowerBound);
    return;
}

class HeldKarp
{
public:
    int size;
    int allVisited;
    int minimum = INT_MAX;
    int* sequence;
    int* path;
    int** edge_matrix;
    int** valueOfSets;
    HeldKarp(int size, int** matrix);
    void startHeldKarp(int mask);
    int D(int mask, int city, int level);
    void displayPath(void);
};

HeldKarp::HeldKarp(int size, int** matrix)
{
    this->size = size;
    allVisited = (1 << size) - 1;
    this->edge_matrix = matrix;
    this->sequence = new int[size];
    this->path = new int[size + 1];
    for (int i = 0; i < size; i++)
        sequence[i] = 0;
    this->valueOfSets = new int* [(1 << size)];
    for (int i = 0; i <= ((1 << size) - 1); i++)
    {
        valueOfSets[i] = new int[size];
    }
    for (int i = 0; i < (1 << size); i++)
    {
        for (int j = 0; j < size; j++)
            valueOfSets[i][j] = -1;
    }
}
void HeldKarp::displayPath()
{
    for (int i = 0; i < size; i++)
        cout << path[i] << " -> ";
    cout << path[size] << endl;
}


int HeldKarp::D(int mask, int prev, int level)
{
    if (mask == allVisited)
    {
        return edge_matrix[prev][0];
    }
    if (valueOfSets[mask][prev] != -1)
    {
        return valueOfSets[mask][prev];
    }
    int actLen = INT_MAX;
    for (int i = 0; i < size; i++)
    {
        if ((mask & (1 << i)) == 0)
        {
            sequence[level] = i;
            int len = edge_matrix[prev][i] + D(mask | (1 << i), i, level + 1);
            actLen = min(len, actLen);
            if (level == 1 && actLen < minimum)
            {
                minimum = actLen;
                for (int i = 0; i < size; i++)
                    path[i] = sequence[i];
                path[size] = sequence[0];
            }
            actLen = min(len, actLen);

        }
    }
    return valueOfSets[mask][prev] = actLen;
}
void HeldKarp::startHeldKarp(int mask)
{
    sequence[0] = 0;
    int length = D(1, 0, 1);
    cout << endl << length << "  " << 100.0 * length / optimum << "%" << endl;
}


int** loadingFile(int* num_of_cities, int* opt_Ham)
{
    string file_name;
    string fileName;
    int size;
    cout << "Wprowadz nazwe pliku:\n";
    cin >> fileName;
    ifstream load(fileName);
    while (load.is_open())
    {
        load >> file_name;
        load >> size;
        *num_of_cities = size;
        int** edge_matrix = new int* [*num_of_cities];
        for (int i = 0; i < *num_of_cities; i++)
        {
            edge_matrix[i] = new int[*num_of_cities];
        }
        if (file_name != fileName)
        {
            cout << "Blad wczytywania" << endl;
        }
        else
        {
            for (int i = 0; i < *num_of_cities; i++)
            {
                for (int j = 0; j < *num_of_cities; j++)
                {
                    load >> edge_matrix[i][j];
                }

            }
            load >> *opt_Ham;
            load.close();
            return edge_matrix;
        }
    }
    return 0;
}

void displayMatrix(int** edge_matrix, int num_of_cities)
{
    for (int i = 0; i < num_of_cities; i++)
    {
        for (int j = 0; j < num_of_cities; j++)
        {
            cout << edge_matrix[i][j] << " ";
        }
        cout << endl;
    }
}


int hamLen(int** edge_matrix, int* num_of_cities)
{
    int len = 0;
    for (int i = 1; i < *num_of_cities; i++)
    {
        len = len + edge_matrix[i - 1][i];
    }
    len = len + edge_matrix[*num_of_cities - 1][0];
    return len;
}


int main()
{
    srand(time(NULL));
    int num_of_cities = 0;
    int opt_Ham = 0;
    int opt = 0;
    int* wsk_num_of_cities = &num_of_cities;
    string file_name;
    int** edge_matrix;
    int* wsk_opt_Ham = &opt_Ham;
    int length;
    cout << "---------------------------------------------------------------------------------" << endl;
    cout << "Wczytywanie instancji\n";
    edge_matrix = loadingFile(wsk_num_of_cities, wsk_opt_Ham);
    optimum = opt_Ham;
    int* sequence = new int[num_of_cities];
    while (1)
    {
        cout << endl << "---------------------------------------------------------------------------------" << endl;
        cout << "Wybierz opcje:\n";
        cout << "1. Wyswietl macierz krawedzi\n";
        cout << "2. Algorytm genetyczny\n";
        cout << "3. Algorytm zachlanny\n";
        cout << "4. Algorytm Brute Force\n";
        cout << "5. Algorytm B&B\n";
        cout << "6. Algorytm Helda Karpa\n";
        cout << "7. Algorytm Symulowanego wyzarzania\n";
        cout << "8. Algorytm Przeszukiwania z zakazami\n";
        cout << "9. Wyswietl dlugosc optymalna cyklu\n";
        cout << "0. Zakoncz\n";
        cin >> opt;
        switch (opt)
        {
        case 0:
        {
            return 0;
            break;
        }
        case 1:
        {
            system("cls");
            displayMatrix(edge_matrix, num_of_cities);
            break;
        }
        case 2:
        {
            system("cls");
            GenethicAlgorith* ga = new GenethicAlgorith(num_of_cities, edge_matrix, num_of_cities * 10, 0.85, 0.1);// num_of_cities * 10
            auto start = chrono::system_clock::now();
            ga->startGA(num_of_cities*10);

            cout << endl << "Dlugosc trasy = " << ga->best.length << " PRD = " << 100.0 * (float)ga->best.length / (float)optimum << "%" << endl;

            ga->displayPath(ga->best.chromosome);
            auto end = chrono::system_clock::now();
            auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
            if (elapsed.count() >= 1000)
            {
                cout << "Czas: " << (float)elapsed.count() / 1000.0 << " s" << endl;
            }
            else
                cout << "Czas: " << elapsed.count() << " ms" << endl;
            break;
        }
        case 3:
        {
            int len = 0;
            system("cls");
            auto start = chrono::system_clock::now();
            sequence = nearest_neigh(edge_matrix, num_of_cities);
            auto end = chrono::system_clock::now();
            auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
            for (int i = 0; i < num_of_cities; i++)
            {
                cout << sequence[i] << " -> ";
                len += edge_matrix[sequence[i]][sequence[(i + 1) % num_of_cities]];
            }
            cout << sequence[0] << endl;
            cout << endl << "Dlugosc trasy = " << len << " PRD = " << 100.0 * (float)len / (float)optimum << "%" << endl;
            cout << "Czas: " << elapsed.count() << " ms" << endl;
            break;
        }
        case 4:
        {
            system("cls");
            auto start = chrono::system_clock::now();
            startBruteF(edge_matrix, wsk_num_of_cities, sequence);
            auto end = chrono::system_clock::now();
            auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
            if (elapsed.count() >= 1000)
            {
                cout << "Czas: " << (float)elapsed.count() / 1000.0 << " s" << endl;
            }
            else
                cout << "Czas: " << elapsed.count() << " ms" << endl;
            for (int i = 0; i < num_of_cities; i++)
                cout << sequence[i] << " -> ";
            cout << sequence[0];
            break;
        }
        case 5:
        {
            system("cls");
            BandBClass* bandb = new BandBClass(num_of_cities, edge_matrix);
            auto start = chrono::system_clock::now();
            bandb->b_and_b_start();
            auto end = chrono::system_clock::now();
            auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
            if (elapsed.count() >= 1000)
            {
                cout << "Czas: " << (float)elapsed.count() / 1000.0 << " s" << endl;
            }
            else
                cout << "Czas: " << elapsed.count() << " ms" << endl;
            bandb->displayPath();
            break;
        }
        case 6:
        {
            system("cls");
            HeldKarp* hk = new HeldKarp(num_of_cities, edge_matrix);
            auto start = chrono::system_clock::now();
            hk->startHeldKarp(0);
            auto end = chrono::system_clock::now();
            auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
            if (elapsed.count() >= 1000)
            {
                cout << "Czas: " << (float)elapsed.count() / 1000.0 << " s" << endl;
            }
            else
                cout << "Czas: " << elapsed.count() << " ms" << endl;
            cout << "Pamiec: " << (1 << num_of_cities) * num_of_cities * sizeof(int) << endl;
            hk->displayPath();
            break;
        }
        case 7:
        {
            system("cls");
            SimulatedAnnealing* sa = new SimulatedAnnealing(num_of_cities, edge_matrix);
            auto start = chrono::system_clock::now();
            sa->startSA(1000000.0, 0.999, 1.0);
            auto end = chrono::system_clock::now();
            auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
            if (elapsed.count() >= 1000)
            {
                cout << "Czas: " << (float)elapsed.count() / 1000.0 << " s" << endl;
            }
            else
                cout << "Czas: " << elapsed.count() << " ms" << endl;
            displayPath(sa->final_path, num_of_cities);
            delete sa;

            break;
        }
        case 8:
        {
            system("cls");
            TabuSearch* ts = new TabuSearch(num_of_cities, edge_matrix);
            auto start = chrono::system_clock::now();
            ts->startTS(10);
            auto end = chrono::system_clock::now();
            auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
            if (elapsed.count() >= 1000)
            {
                cout << "Czas: " << (float)elapsed.count() / 1000.0 << " s" << endl;
            }
            else
                cout << "Czas: " << elapsed.count() << " ms" << endl;
            displayPath(ts->final_path, num_of_cities);
            break;
        }
        case 9:
        {
            system("cls");
            cout << "Dlugosc optymalna cyklu = " << opt_Ham << endl;
            break;
        }
        default:
            break;
        }
    }
}
