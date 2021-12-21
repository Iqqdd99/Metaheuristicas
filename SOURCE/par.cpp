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
    for(int i = 0; i < clusters.size(); i++)
    {
        if(clusters[i] == cluster_item){
            count++;
        }
    }
    
    
    return count == 1;
}

float heuristicFunction(const vector<int> &clusters, const vector<vector<float>> &data, const vector<vector<int>> &constrictions, int number_clusters, float lambda)
{
    return generalPartitionDesviation(clusters, data, number_clusters) + (lambda * infeasability(clusters, constrictions));
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
    int number_clusters = stoi(argv[3]);
    int count = 0;
    int seed = stoi(argv[4]);
    int change_count = 0;
    float lambda;

    Set_random(seed);


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

    cout << " ===== Semilla: " << seed << " Data: " << argv[1] << " Constrictions: " << argv[2] << " =====" << endl;

    lambda = generateLambdaValue(data, constrictions);

    float stop_value = data.size() * 0.025;
    cout << "\tBúsqueda Local" << endl;
    auto start = chrono::steady_clock::now();
    for (int i = 0; i < data.size(); i++)
    {
        solution.push_back(Randint(0, number_clusters - 1));
    }

    vector<int> index;
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
                solution[index[i]] = best_cluster;
            }
        }

        change_count = 0;
        for (int i = 0; i < solution.size(); i++)
        {
            if (old_solution[i] != solution[i])
            {
                change_count++;
            }
        }
    } while (stop_value < change_count); //change_count > stop_value);
    auto end = chrono::steady_clock::now();

    printSolution(solution, data);
    cout << "Tasa_C " << generalPartitionDesviation(solution, data, number_clusters) << endl;
    cout << "Tasa_inf " << infeasability(solution, constrictions) << endl;
    cout << "Agregado " << heuristicFunction(solution, data, constrictions, number_clusters, lambda) << endl;
    cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

    solution.clear();
    index.clear();

    cout << "\tGreedy" << endl;
    start = chrono::steady_clock::now();

    vector<vector<float>> centroids(number_clusters, vector<float>(data[0].size(), 0));
    int feasable_clussters = number_clusters / 2;
    vector<pair<int, int>> sorted_infes;

    for (int i = 0; i < data.size(); i++)

    {
        solution.push_back(Randint(0, number_clusters - 1));
        }

        for (int i = 0; i < solution.size(); i++)
        {
            index.push_back(i);
        }

        for (int i = 0; i < data[0].size(); i++)
        {
            float higher_value = 0, lower_value = data.at(0).at(0);
            for (int j = 0; j < data.size(); j++)
            {
                if (data[j][i] > higher_value)
                {
                    higher_value = data[j][i];
                }

                if (data[j][i] < lower_value)
                {
                    lower_value = data[j][i];
                }
            }
            for (int j = 0; j < centroids.size(); j++)
            {
                centroids[j][i] = Randfloat(lower_value, higher_value);
            }

        }

        random_shuffle(index.begin(), index.end());
        do
        {
            for (int i = 0; i < index.size(); i++)
            {
                if (!lastElement(solution, index[i]))
                {
                    for (int j = 0; j < number_clusters; j++)
                    {
                        vector<int> partial_solution = solution;
                        partial_solution[index[i]] = j;
                        sorted_infes.push_back( make_pair(infeasability(partial_solution, constrictions),j) );
                    }
                    sort(sorted_infes.begin(),sorted_infes.end());

                    float nearest_distance = EuclideanDistance(centroids[0],data[index[i]]);
                    int nearest_centroid = 0;

                    for(int j = 1; j < number_clusters-feasable_clussters; j++)
                    {
                        float distance = EuclideanDistance(centroids[sorted_infes[j].second],data[index[i]]);
                        if(distance < nearest_distance){
                            nearest_distance = distance;
                            nearest_centroid = sorted_infes[j].second;
                        }
                    }
                    
                    solution[index[i]] = nearest_centroid;
                    sorted_infes.clear();
                }
            }
            
            for (int j = 0; j < number_clusters; j++)
            {
                centroids[j] = clusterCentroid(solution, data, j);
            }
            count++;
        } while (count <500);
        end = chrono::steady_clock::now();

        printSolution(solution, data);
        cout << "Tasa_C " << generalPartitionDesviation(solution, data, number_clusters) << endl;
        cout << "Tasa_inf " << infeasability(solution, constrictions) << endl;
        cout << "Agregado " << heuristicFunction(solution, data, constrictions, number_clusters, lambda) << endl;
        cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

        solution.clear();
        index.clear();

    return 0;

}
