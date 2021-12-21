//Compilacion: g++  -std=c++11 -O3  interaccion_grav.cpp random.cpp -o grav
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <stdlib.h>
#include <vector>
#include <utility>
#include <chrono>
#include <unistd.h>
#include <algorithm>
#include <math.h>
#include "random.h"

using namespace std;

float EuclideanDistance(const vector<float> &point_a, const vector<float> &point_b)
{
    float distance = 0;
    for (int i = 0; i < point_a.size(); i++)
    {
        float aux = point_a.at(i) - point_b.at(i);
        distance += aux * aux;
    }
    return sqrt(distance);
}

vector<float> clusterCentroid(const vector<int> &clusters, const vector<vector<float>> &data, int cluster)
{
    vector<float> cent(data[0].size(), 0);
    int size = 0;
    for (int i = 0; i < clusters.size(); i++)
    {
        if (clusters[i] == cluster)
        {
            size++;
            for (int j = 0; j < data[0].size(); j++)
            {
                cent.at(j) += data[i][j];
            }
        }
    }

    for (int i = 0; i < data[0].size(); i++)
    {
        if (cent[i] != 0)
        {
            cent[i] /= size;
        }
    }

    return cent;
}

float calculateIntraClusterDistance(const vector<int> &clusters, const vector<vector<float>> &data, int cluster)
{
    vector<float> centroid = clusterCentroid(clusters, data, cluster);
    float Distance = 0;
    int size = 0;
    for (int i = 0; i < clusters.size(); i++)
    {
        if (clusters[i] == cluster)
        {
            size++;
            Distance += EuclideanDistance(data[i], centroid);
        }
    }

    return Distance / size;
}
float generateLambdaValue(const vector<vector<float>> &data, const vector<vector<int>> &constrictions)
{

    int constr = 0;
    float max_distance = 0;
    for (int i = 0; i < constrictions.size(); i++)
    {
        for (int j = i + 1; j < constrictions[0].size(); j++)
        {
            if (constrictions[i][j] != 0)
            {
                constr++;
            }
        }
    }

    for (int i = 0; i < data.size(); i++)
    {
        for (int j = i + 1; j < data.size(); j++)
        {
            float distance = EuclideanDistance(data[i], data[j]);
            if (distance > max_distance)
            {
                max_distance = distance;
            }
        }
    }

    return max_distance / constr;
}

float generalPartitionDesviation(const vector<int> &clusters, const vector<vector<float>> &data, int number_clusters)
{
    float desviation = 0;
    for (int i = 0; i < number_clusters; i++)
    {
        desviation += calculateIntraClusterDistance(clusters, data, i);
    }
    return desviation / number_clusters;
}

int infeasability(const vector<int> &clusters, const vector<vector<int>> &constrictions)
{
    int value = 0;

    for (int i = 0; i < constrictions.size(); i++)
    {
        for (int j = i + 1; j < constrictions[0].size(); j++)
        {
            if (constrictions[i][j] != 0)
            {
                if (constrictions[i][j] == 1)
                {
                    if (clusters[j] != clusters[i])
                    {
                        value++;
                    }
                }
                else
                {
                    if (clusters[j] == clusters[i])
                    {
                        value++;
                    }
                }
            }
        }
    }

    return value;
}

bool lastElement(const vector<int> &clusters, int item)
{
    int cluster_item = clusters[item];
    int count = 0;
    for (int i = 0; i < clusters.size(); i++)
    {
        if (clusters[i] == cluster_item)
        {
            count++;
        }
    }

    return count == 1;
}

float heuristicFunction(const vector<int> &clusters, const vector<vector<float>> &data, const vector<vector<int>> &constrictions, int number_clusters, float lambda)
{
    //float gen = generalPartitionDesviation(clusters, data, number_clusters);
    //float infes =  infeasability(clusters, constrictions);
    return generalPartitionDesviation(clusters, data, number_clusters) + (lambda * infeasability(clusters, constrictions));
    //return gen + (lambda * infes);
}

void printSolution(const vector<int> &clusters, const vector<vector<float>> &data)
{
    for (int i = 0; i < clusters.size(); i++)
    {
        cout << clusters[i] << " ";
        for (int j = 0; j < data[0].size(); j++)
        {
            cout << data[i][j] << " ";
        }
        cout << endl;
    }
}

int correctSolution(const vector<int> &clusters, int number_clusters)
{
    vector<bool> check(number_clusters, false);
    vector<bool> correct_solution(number_clusters, true);
    for (int i = 0; i < clusters.size(); i++)
    {
        check[clusters[i]] = true;
    }

    if (check == correct_solution)
    {
        return -1;
    }
    else
    {
        for (int i = 0; i < check.size(); i++)
        {
            if (!check[i])
                return i;
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        cout << "error in the arguments. Correct usage: ./par file.dat file_const.dat number_clusters seed_value" << endl;
        return -1;
    }

    vector<vector<float>> data;
    vector<vector<int>> constrictions;
    vector<int> solution;
    vector<vector<int>> solutions;
    vector<float> solutions_values;
    vector<int> bestSolution;
    vector<int> index;
    float bestValue;
    int indexBestValue;
    int number_clusters = stoi(argv[3]);
    int nObjetivos = 10;
    int seed = stoi(argv[4]);
    int heuristicEvaluations;
    int maxHeuristicEvaluations = 10000;
    float lambda;
    float G = 6.674e-11;
    chrono::steady_clock::time_point start;
    chrono::steady_clock::time_point end;

    Set_random(seed);

    /*
     *  Lectura de datos y restricciones   
     **/
    ifstream File;
    File.open(argv[1]);

    if (File.is_open())
    {
        string input;
        while (!File.eof())
        {
            getline(File, input);
            stringstream ss(input);
            vector<float> point;
            while (ss.good())
            {
                string value;
                getline(ss, value, ',');
                point.push_back(stof(value));
            }

            data.push_back(point);
        }
    }
    File.close();

    File.open(argv[2]);
    if (File.is_open())
    {
        string input;
        while (!File.eof())
        {
            getline(File, input);
            stringstream ss(input);
            vector<int> point;
            while (ss.good())
            {
                string value;
                getline(ss, value, ',');
                if (value != "")
                {
                    point.push_back(stoi(value));
                }
            }
            constrictions.push_back(point);
        }
    }
    File.close();

    /*
     *  Generación de valor Lambda para la función de evaluación
     * */

    lambda = generateLambdaValue(data, constrictions);

    cout << " ===== Semilla: " << seed << " Data: " << argv[1] << " Constrictions: " << argv[2] << " =====" << endl;

    cout << "\tBúsqueda Mediante gravitación" << endl;

    start = chrono::steady_clock::now();

    heuristicEvaluations = 0;

    /*
     * Generación de solución inicial
     **/

    do
    {
        solution.clear();
        for (int i = 0; i < data.size(); i++)
        {
            solution.push_back(Randint(0, number_clusters - 1));
        }
    } while (correctSolution(solution, number_clusters) != -1);

     /*
     * Generación de índice
     **/

    for (int i = 0; i < solution.size(); i++)
    {
        index.push_back(i);
    }

     /*
     * Bucle principal del algoritmo
     **/
    do
    {
        random_shuffle(index.begin(), index.end());
        for (int i = 0; i < index.size(); i++)
        {
            /*
             * Calculo de las masas de los clusters
             **/
            vector<int> mass(number_clusters);
            for (int j = 0; j < index.size(); j++)
            {
                mass[solution[index[j]]]++;
            }

            /*
             * Calculo de los centroides de los clusters
             **/
            vector<vector<float>> centroids(number_clusters);
            for (int j = 0; j < number_clusters; j++)
            {
                centroids[j] = clusterCentroid(solution, data, j);
            }

            if (!lastElement(solution, index[i]))
            {
                /*
                 * Calculo de la fuerza de interacción gravitatoria del punto con los clusters
                 **/
                float max_force = G * (mass[0] / EuclideanDistance(centroids[0], data[index[i]]));
                int best_cluster = 0;

                for (int j = 1; j < number_clusters; j++)
                {
                    float point_force = G * (mass[j] / EuclideanDistance(centroids[j], data[index[i]]));
                    if (point_force > max_force)
                    {
                        max_force = point_force;
                        best_cluster = j;
                    }
                }
                heuristicEvaluations++;
                solution[index[i]] = best_cluster;
            }
        }

    } while (heuristicEvaluations < maxHeuristicEvaluations);
    
    end = chrono::steady_clock::now();

    cout << generalPartitionDesviation(solution, data, number_clusters)
         << "," << infeasability(solution, constrictions)
         << "," << heuristicFunction(solution, data, constrictions, number_clusters, lambda)
         << "," << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;

    index.clear();
    solution.clear();

    cout << "\tBúsqueda Mediante Gravitación multiarranque" << endl;
    start = chrono::steady_clock::now();

    heuristicEvaluations = 0;

    /*
     * Generación de índice
     **/
    for (int i = 0; i < data.size(); i++)
    {
        index.push_back(i);
    }

    /*
     * Bucle que genera los diferentes arranques
     **/
    for (int solucion_i = 0; solucion_i < nObjetivos; solucion_i++)
    {
        /*
         * Generación de solución inicial
         **/
        do
        {
            solution.clear();
            for (int i = 0; i < data.size(); i++)
            {
                solution.push_back(Randint(0, number_clusters - 1));
            }
        } while (correctSolution(solution, number_clusters) != -1);

        /*
         * Bucle principal del algoritmo
         **/
        do
        {
            random_shuffle(index.begin(), index.end());
            for (int i = 0; i < index.size(); i++)
            {
                /*
                 * Calculo de las masas de los clusters
                 **/
                vector<int> mass(number_clusters);
                for (int j = 0; j < index.size(); j++)
                {
                    mass[solution[index[j]]]++;
                }

                /*
                 * Calculo de los centroides de los clusters
                 **/
                vector<vector<float>> centroids(number_clusters);
                for (int j = 0; j < number_clusters; j++)
                {
                    centroids[j] = clusterCentroid(solution, data, j);
                }

                if (!lastElement(solution, index[i]))
                {
                    /*
                     * Calculo de la fuerza de interacción gravitatoria del punto con los clusters
                     **/
                    float max_force = G * (mass[0] / EuclideanDistance(centroids[0], data[index[i]]));
                    int best_cluster = 0;

                    for (int j = 1; j < number_clusters; j++)
                    {
                        float point_force = G * (mass[j] /  EuclideanDistance(centroids[j], data[index[i]]));
                        if (point_force > max_force)
                        {
                            max_force = point_force;
                            best_cluster = j;
                        }
                    }
                                    heuristicEvaluations++;
                    solution[index[i]] = best_cluster;
                }
            }

        } while (heuristicEvaluations < maxHeuristicEvaluations);

        solutions.push_back(solution);
        solutions_values.push_back(heuristicFunction(solution, data, constrictions, number_clusters, lambda));
    }

    bestValue = solutions_values[0];
    indexBestValue = 0;

    /*
     * Selección de la mejor solución
     * */
    for (int solution_i = 1; solution_i < nObjetivos; solution_i++)
    {
        if (solutions_values[solution_i] < bestValue)
        {
            bestValue = solutions_values[solution_i];
            indexBestValue = solution_i;
        }
    }

    solution = solutions[indexBestValue];

    solutions.clear();
    index.clear();
    solutions_values.clear();

    end = chrono::steady_clock::now();

    cout << generalPartitionDesviation(solution, data, number_clusters)
         << "," << infeasability(solution, constrictions)
         << "," << heuristicFunction(solution, data, constrictions, number_clusters, lambda)
         << "," << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;

    solution.clear();
}
