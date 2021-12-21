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

    cout << "\tBúsqueda Local Multiarranque Básica" << endl;

    start = chrono::steady_clock::now();
    for(int solucion_i = 0; solucion_i < nObjetivos; solucion_i++){

        heuristicEvaluations = 0;

        for (int i = 0; i < data.size(); i++)
        {
            solution.push_back(Randint(0, number_clusters - 1));
        }

        for (int i = 0; i < solution.size(); i++)
        {
            index.push_back(i);
        }

        do
        {
            random_shuffle(index.begin(), index.end());
            vector<int> old_solution = solution;
            for (int i = 0; i < index.size(); i++)
            {

                if (!lastElement(solution, index[i]))
                {
                    vector<int> partial_solution = solution;
                    partial_solution[index[i]] = 0;
                    float min_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                    int best_cluster = 0;
                    for (int j = 1; j < number_clusters; j++)
                    {
                        vector<int> partial_solution = solution;
                        partial_solution[index[i]] = j;
                        float function_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                        if (min_value > function_value)
                        {
                            min_value = function_value;
                            best_cluster = j;
                        }
                    }
                    heuristicEvaluations+=number_clusters;
                    solution[index[i]] = best_cluster;
                }
            }

        } while (heuristicEvaluations < maxHeuristicEvaluations);

        solutions.push_back(solution);
        solutions_values.push_back(heuristicFunction(solution, data, constrictions, number_clusters, lambda));

        index.clear();
        solution.clear();

    }

    bestValue = solutions_values[0];
    indexBestValue = 0;

    for(int solucion_i = 1; solucion_i < nObjetivos; solucion_i++){
        if(solutions_values[solucion_i] < bestValue){
            bestValue = solutions_values[solucion_i];
            indexBestValue = solucion_i;
        }
    }

    solution = solutions[indexBestValue];

    solutions.clear();
    solutions_values.clear();

    end = chrono::steady_clock::now();

    // printSolution(best_chromosome, data);
    // cout << "Tasa_C " << generalPartitionDesviation(solution, data, number_clusters) << endl;
    // cout << "Tasa_inf " << infeasability(solution, constrictions) << endl;
    // cout << "Agregado " << heuristicFunction(solution, data, constrictions, number_clusters, lambda) << endl;
    // cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

    cout << generalPartitionDesviation(solution, data, number_clusters)
        << "," << infeasability(solution, constrictions)
        << "," << heuristicFunction(solution, data, constrictions, number_clusters, lambda)
        << "," << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;

    index.clear();
    solution.clear();

    cout << "\tBúsqueda Local Reiterada" << endl;

    start = chrono::steady_clock::now();

    heuristicEvaluations = 0;

    for (int i = 0; i < data.size(); i++)
    {
        solution.push_back(Randint(0, number_clusters - 1));
    }

    for (int i = 0; i < solution.size(); i++)
    {
        index.push_back(i);
    }

    do
    {
        random_shuffle(index.begin(), index.end());
        vector<int> old_solution = solution;
        for (int i = 0; i < index.size(); i++)
        {

            if (!lastElement(solution, index[i]))
            {
                vector<int> partial_solution = solution;
                partial_solution[index[i]] = 0;
                float min_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                int best_cluster = 0;
                for (int j = 1; j < number_clusters; j++)
                {
                    vector<int> partial_solution = solution;
                    partial_solution[index[i]] = j;
                    float function_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                    if (min_value > function_value)
                    {
                        min_value = function_value;
                        best_cluster = j;
                    }
                }
                heuristicEvaluations += number_clusters;
                solution[index[i]] = best_cluster;
            }
        }

    } while (heuristicEvaluations < maxHeuristicEvaluations);

    solutions.push_back(solution);
    solutions_values.push_back(heuristicFunction(solution, data, constrictions, number_clusters, lambda));

    for(int solucion_i = 1; solucion_i < nObjetivos; solucion_i++){

        heuristicEvaluations = 0;

        //todo mutacion
        vector<int> positions(solution.size()*0.1);
        for(int i=0; i<positions.size(); i++)
        {
            positions[i] = Randint(0,solution.size()-1);
        }

        for(int i=0; i<positions.size(); i++){
            solution[positions[i]] = Randint(0,number_clusters-1);
        }

        do
        {
            random_shuffle(index.begin(), index.end());
            vector<int> old_solution = solution;
            for (int i = 0; i < index.size(); i++)
            {

                if (!lastElement(solution, index[i]))
                {
                    vector<int> partial_solution = solution;
                    partial_solution[index[i]] = 0;
                    float min_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                    int best_cluster = 0;
                    for (int j = 1; j < number_clusters; j++)
                    {
                        vector<int> partial_solution = solution;
                        partial_solution[index[i]] = j;
                        float function_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                        if (min_value > function_value)
                        {
                            min_value = function_value;
                            best_cluster = j;
                        }
                    }
                    heuristicEvaluations+=number_clusters;
                    solution[index[i]] = best_cluster;
                }
            }

        } while (heuristicEvaluations < maxHeuristicEvaluations);

        solutions.push_back(solution);
        solutions_values.push_back(heuristicFunction(solution, data, constrictions, number_clusters, lambda));

    }

    bestValue = solutions_values[0];
    indexBestValue = 0;

    for(int solution_i = 1; solution_i < nObjetivos; solution_i++){
        if(solutions_values[solution_i] < bestValue){
            bestValue = solutions_values[solution_i];
            indexBestValue = solution_i;
        }
    }

    solution = solutions[indexBestValue];

    solutions.clear();
    index.clear();
    solutions_values.clear();

    end = chrono::steady_clock::now();

    // printSolution(best_chromosome, data);
    // cout << "Tasa_C " << generalPartitionDesviation(solution, data, number_clusters) << endl;
    // cout << "Tasa_inf " << infeasability(solution, constrictions) << endl;
    // cout << "Agregado " << heuristicFunction(solution, data, constrictions, number_clusters, lambda) << endl;
    // cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

    cout << generalPartitionDesviation(solution, data, number_clusters)
        << "," << infeasability(solution, constrictions)
        << "," << heuristicFunction(solution, data, constrictions, number_clusters, lambda)
        << "," << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;

    solution.clear();

    cout << "\tEnfriamiento simulado" << endl;

    float max_vecinos = 10 * data.size();
    float max_exitos = 0.1 * max_vecinos;
    int exitos;
    int vecinos;
    bool cambios;
    float mValue = 100000 / max_vecinos;
    float beta;
    float initialTemp;
    float finalTemp = 0.01;
    float currentTemp;
    float fi = 0.3;
    float mu = 0.3;

    start = chrono::steady_clock::now();
    heuristicEvaluations = 0;
    maxHeuristicEvaluations = 100000;

    for (int i = 0; i < data.size(); i++)
    {
        solution.push_back(Randint(0, number_clusters - 1));
    }

    for (int i = 0; i < solution.size(); i++)
    {
        index.push_back(i);
    }

    bestSolution = solution;
    bestValue = generalPartitionDesviation(bestSolution, data, number_clusters);
    currentTemp = initialTemp = (mu * bestValue) / (-log(fi));
    beta = (initialTemp-finalTemp)/(mValue*finalTemp*initialTemp);
    cambios = false;

    do
    {
        exitos = 0;
        vecinos = 0;
        int i = 0;
        int old_cluster;
        random_shuffle(index.begin(), index.end());
        vector<int> old_solution = solution;
        do
        {
            if (!lastElement(solution, index[i]))
            {
                old_cluster = solution[index[i]];
                vector<int> partial_solution = solution;
                partial_solution[index[i]] = 0;
                float min_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                int best_cluster = 0;
                for (int j = 1; j < number_clusters; j++)
                {
                    vector<int> partial_solution = solution;
                    partial_solution[index[i]] = j;
                    float function_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                    if (min_value > function_value)
                    {
                        min_value = function_value;
                        best_cluster = j;
                    }
                }
                solution[index[i]] = best_cluster;

                if (best_cluster != old_cluster)
                {
                    cambios = true;
                    exitos++;
                }

                if (min_value < bestValue)
                {
                    bestValue = min_value;
                    bestSolution = solution;
                }
            }
            else
                cambios = true;

            i++;
            vecinos++;
            if (i == index.size())
            {
                i = 0;
               // random_shuffle(index.begin(), index.end());
            }
        } while (exitos < max_exitos && vecinos < max_vecinos);
        currentTemp = (currentTemp) / (1 + (beta * currentTemp));
        // currentTemp = 0.9 * currentTemp;
    } while (currentTemp > finalTemp && cambios);

    end = chrono::steady_clock::now();

    // printSolution(best_chromosome, data);
    // cout << "Tasa_C " << generalPartitionDesviation(solution, data, number_clusters) << endl;
    // cout << "Tasa_inf " << infeasability(solution, constrictions) << endl;
    // cout << "Agregado " << heuristicFunction(solution, data, constrictions, number_clusters, lambda) << endl;
    // cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

    cout << generalPartitionDesviation(bestSolution, data, number_clusters)
         << "," << infeasability(bestSolution, constrictions)
         << "," << heuristicFunction(bestSolution, data, constrictions, number_clusters, lambda)
         << "," << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;

    solution.clear();
    index.clear();
    bestSolution.clear();

    cout << "\tEnfriamiento Simulado Reiterado " << endl;

    start = chrono::steady_clock::now();

    maxHeuristicEvaluations = 10000;
    heuristicEvaluations = 0;
    mValue = 10000 / max_vecinos;

    for (int i = 0; i < data.size(); i++)
    {
        solution.push_back(Randint(0, number_clusters - 1));
    }

    for (int i = 0; i < solution.size(); i++)
    {
        index.push_back(i);
    }

    bestSolution = solution;
    bestValue = generalPartitionDesviation(bestSolution, data, number_clusters);
    currentTemp = initialTemp = (mu * bestValue) / (-log(fi));

    beta = (initialTemp - finalTemp) / (mValue * initialTemp * finalTemp);
    cambios = false;

    do
    {
        exitos = 0;
        vecinos = 0;
        int i = 0;
        int old_cluster;
        random_shuffle(index.begin(), index.end());
        vector<int> old_solution = solution;
        do
        {
            if (!lastElement(solution, index[i]))
            {
                old_cluster = solution[index[i]];
                vector<int> partial_solution = solution;
                partial_solution[index[i]] = 0;
                float min_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                int best_cluster = 0;
                for (int j = 1; j < number_clusters; j++)
                {
                    vector<int> partial_solution = solution;
                    partial_solution[index[i]] = j;
                    float function_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                    if (min_value > function_value)
                    {
                        min_value = function_value;
                        best_cluster = j;
                    }
                }
                solution[index[i]] = best_cluster;

                if (best_cluster != old_cluster)
                {
                    cambios = true;
                    exitos++;
                }

                if (min_value < bestValue)
                {
                    bestValue = min_value;
                    bestSolution = solution;
                }
            }
            else
                cambios = true;

            i++;
            vecinos++;
            if (i == index.size())
            {
                i = 0;
            }
        } while (exitos < max_exitos && vecinos < max_vecinos);
        currentTemp = (currentTemp) / (1 + (beta * currentTemp));
        // currentTemp = 0.9 * currentTemp;
    } while (currentTemp > finalTemp && cambios);

    solutions.push_back(bestSolution);
    solutions_values.push_back(heuristicFunction(bestSolution, data, constrictions, number_clusters, lambda));  

    solution = bestSolution;

    for (int solucion_i = 1; solucion_i < nObjetivos; solucion_i++)
    {

        heuristicEvaluations = 0;

        //todo mutacion
            vector<int> positions(solution.size() * 0.1);
            for (int i = 0; i < positions.size(); i++)
            {
                positions[i] = Randint(0, solution.size() - 1);
            }

            int aux;
            for (int i = 0; i < positions.size(); i++)
            {
                solution[positions[i]] = Randint(0, number_clusters - 1);
            }

            bestSolution = solution;
            bestValue = generalPartitionDesviation(bestSolution, data, number_clusters);
            currentTemp = (mu * bestValue) / (-log(fi));

            currentTemp = initialTemp = (mu * bestValue) / (-log(fi));

            beta = (initialTemp - finalTemp) / (mValue * initialTemp * finalTemp);
            cambios = false;

            do
            {
                exitos = 0;
                vecinos = 0;
                int i = 0;
                int old_cluster;
                random_shuffle(index.begin(), index.end());
                vector<int> old_solution = solution;
                do
                {
                    if (!lastElement(solution, index[i]))
                    {
                        old_cluster = solution[index[i]];
                        vector<int> partial_solution = solution;
                        partial_solution[index[i]] = 0;
                        float min_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                        int best_cluster = 0;
                        for (int j = 1; j < number_clusters; j++)
                        {
                            vector<int> partial_solution = solution;
                            partial_solution[index[i]] = j;
                            float function_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                            if (min_value > function_value)
                            {
                                min_value = function_value;
                                best_cluster = j;
                            }
                        }
                        solution[index[i]] = best_cluster;

                        if (best_cluster != old_cluster)
                        {
                            cambios = true;
                            exitos++;
                        }

                        if (min_value < bestValue)
                        {
                            bestValue = min_value;
                            bestSolution = solution;
                        }
                    }
                    else
                        cambios = true;

                    i++;
                    vecinos++;
                    if (i == index.size())
                    {
                        i = 0;
                      //  random_shuffle(index.begin(), index.end());
                    }
                     
                } while (exitos < max_exitos && vecinos < max_vecinos);
                currentTemp = (currentTemp) / (1 + (beta * currentTemp));
                // currentTemp = 0.9 * currentTemp;
            } while (currentTemp > finalTemp && cambios);

            solutions.push_back(bestSolution);
            solutions_values.push_back(heuristicFunction(bestSolution, data, constrictions, number_clusters, lambda));
            solution = bestSolution;
    }

    bestValue = solutions_values[0];
    indexBestValue = 0;

    for (int solution_i = 1; solution_i < nObjetivos; solution_i++)
    {
        if (solutions_values[solution_i] < bestValue)
        {
            bestValue = solutions_values[solution_i];
            indexBestValue = solution_i;
        }
    }

    solution = solutions[indexBestValue];

    index.clear();

    end = chrono::steady_clock::now();

    // printSolution(best_chromosome, data);
    // cout << "Tasa_C " << generalPartitionDesviation(solution, data, number_clusters) << endl;
    // cout << "Tasa_inf " << infeasability(solution, constrictions) << endl;
    // cout << "Agregado " << heuristicFunction(solution, data, constrictions, number_clusters, lambda) << endl;
    // cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

    cout << generalPartitionDesviation(solution, data, number_clusters)
         << "," << infeasability(solution, constrictions)
         << "," << heuristicFunction(solution, data, constrictions, number_clusters, lambda)
         << "," << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;

    solution.clear();
    }
