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

vector<vector<int>> generateIntitialChromosomes(const vector<vector<float>> &data, int number_clusters, int number_chromosomes)
{
    vector<vector<int>> chromosomes(number_chromosomes, vector<int>(data.size(), -1));
    for (int i = 0; i < chromosomes.size(); i++)
    {
        do
        {
            for (int j = 0; j < chromosomes[0].size(); j++)
            {
                chromosomes[i][j] = Randint(0, number_clusters - 1);
            }
        } while (correctSolution(chromosomes[i], number_clusters) != -1);
    }
    return chromosomes;
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
    return generalPartitionDesviation(clusters, data, number_clusters) + (lambda * infeasability(clusters, constrictions));
}

void fixOperator(vector<int> &clusters, int number_clusters){
    
}

void mutationOperator(vector<vector<int>> &chromosomes, int number_clusters){
    for(int i=0; i<chromosomes.size(); i++){
        if(Randfloat(0,1) <= 0.001){
            int pos = Randint(0,chromosomes[0].size()-1);
            int new_cluster = Randint(0,number_clusters-1);
            while(new_cluster == chromosomes[i][pos]){
                 new_cluster = Randint(0,number_clusters-1);
            }
            chromosomes[i][pos] = new_cluster;
        }
    }
}

void uniformCrossoverOperator(vector<vector<int>> &chromosomes, int number_clusters, float prob)
{
    vector<int> positions(chromosomes[0].size()/2);
    int aux;
    for(int i=0; i<positions.size(); i++)
    {
        positions[i] = Randint(0,chromosomes[0].size()-1);
    }

    for(int i=0; i<chromosomes.size(); i+=2){
        if(Randfloat(0,1) <= prob){
            for(int j=0; j<positions.size(); j++){
                aux = chromosomes[i+1][j];
                chromosomes[i+1][j] = chromosomes[i][j];
                chromosomes[i][j] = aux;
            }
        }
    }
    
}

void FixedSegmentCrossoverOperator(vector<vector<int>> &chromosomes, int number_clusters, float prob)
{
    int pos_1 = Randint(0,chromosomes[0].size()-1);
    int pos_2 = Randint(0,chromosomes[0].size()-1);
   
    if(pos_1 > pos_2){
        int aux = pos_2;
        pos_2 = pos_1;
        pos_1 = aux;
    }

    int seg = (pos_2-pos_1)/2;
    vector<int> positions(seg);

    for(int i=0; i<positions.size(); i++)
    {
        positions[i] = Randint(pos_1,pos_2-1);
    }

    for(int i=0; i<chromosomes.size(); i+=2){
        if(Randfloat(0,1) <= prob){
            for(int j=0; j<pos_1; j++){
                chromosomes[i+1][j] = chromosomes[i][j];
            }

            for(int j=0; j<seg; j++){
                chromosomes[i+1][positions[j]] = chromosomes[i][positions[j]];
            }

            for(int j=pos_2; j<chromosomes[0].size(); j++){
                chromosomes[i+1][j] = chromosomes[i][j];
            }
        }
    }
    
 
}

/* 
 * Seleccionará una población de padres del mismo tamaño que la población genética usando un torneo binario con contendientes aleatorios
 * Se realiza el cruce y se cambia la antigua población por la nueva manteniendo el mejor cromosoma de la antigua y eliminando el peor de la nueva.
 * */
vector<vector<int>> UniformElitistEvolutionScheme(const vector<vector<int>> &chromosomes, const vector<vector<float>> &data, const vector<vector<int>> &constrictions, int number_clusters, float lambda, int &heuristicEvaluations)
{
    vector<vector<int>> selected_chromosomes(chromosomes.size());

    //calculo del valor heurístico del primer valor
    vector<int> best_chromosome = chromosomes[0];
    float best_chromosome_value = heuristicFunction(chromosomes[0],data,constrictions,number_clusters,lambda);
    
    int first_wannabe;
    int second_wannabe;
    float prob = 0.7;

    //selección de los primeros contendientes
    do{
        first_wannabe = Randint(0, chromosomes.size() - 1);
        second_wannabe = Randint(0, chromosomes.size() - 1);
    } while(first_wannabe == second_wannabe);

    float first_wannabe_value = heuristicFunction(chromosomes[first_wannabe],data,constrictions,number_clusters,lambda);
    float second_wannabe_value = heuristicFunction(chromosomes[second_wannabe],data,constrictions,number_clusters,lambda);


    if(first_wannabe_value < second_wannabe_value)
    {
        selected_chromosomes[0] = chromosomes[first_wannabe];
    } else {
        selected_chromosomes[0] = chromosomes[second_wannabe];
    }

    // En el mismo bucle calculamos el mejor valor de los chromosomas precursores y hacemos los torneos binarios
    for(int i=1; i<chromosomes.size(); i++)
    {
        //torneo binario
        do
        {
            first_wannabe = Randint(0, chromosomes.size() - 1);
            second_wannabe = Randint(0, chromosomes.size() - 1);
        } while (first_wannabe == second_wannabe);

        first_wannabe_value = heuristicFunction(chromosomes[first_wannabe],data,constrictions,number_clusters,lambda);
        second_wannabe_value = heuristicFunction(chromosomes[second_wannabe],data,constrictions,number_clusters,lambda);

        if(first_wannabe_value < second_wannabe_value)
        {
            selected_chromosomes[i] = chromosomes[first_wannabe];
        }
        else
        {
            selected_chromosomes[i] = chromosomes[second_wannabe];
        }

        //busqueda del mejor valor
        float heuristic_value = heuristicFunction(chromosomes[i],data,constrictions,number_clusters,lambda);
        if(heuristic_value < best_chromosome_value){
            best_chromosome = chromosomes[i];
            best_chromosome_value = heuristic_value;
        }
    }

    uniformCrossoverOperator(selected_chromosomes, number_clusters, prob);
    //selected_chromosomes[worst_chromosome_index] = best_chromosome;

    //Aplicamos mutacion
    mutationOperator(selected_chromosomes, number_clusters);

    //Operador de arreglo
    for(int i = 0; i<selected_chromosomes.size(); i++){
        int isCorrect = correctSolution(selected_chromosomes[i],number_clusters);
        while(isCorrect != -1){
            selected_chromosomes[i][Randint(0,selected_chromosomes[i].size()-1)] = isCorrect;
            isCorrect = correctSolution(selected_chromosomes[i],number_clusters);
        }
    }


    int worst_position = 0;
    float worst_chromosome_value = heuristicFunction(selected_chromosomes[0],data,constrictions,number_clusters,lambda);
 
    for(int i = 1; i<selected_chromosomes.size(); i++){
        float heuristic_value = heuristicFunction(selected_chromosomes[i],data,constrictions,number_clusters,lambda);
        if(heuristic_value > worst_chromosome_value){
            worst_position = i;
            worst_chromosome_value = heuristic_value;
        }
    }

    selected_chromosomes[worst_position] = best_chromosome;

    heuristicEvaluations += chromosomes.size()*4; 

    return selected_chromosomes;
}

vector<vector<int>> FixedSegmentElitistEvolutionScheme(const vector<vector<int>> &chromosomes, const vector<vector<float>> &data, const vector<vector<int>> &constrictions, int number_clusters, float lambda, int &heuristicEvaluations)
{
    vector<vector<int>> selected_chromosomes(chromosomes.size());

    //calculo del valor heurístico del primer valor
    vector<int> best_chromosome = chromosomes[0];
    float best_chromosome_value = heuristicFunction(chromosomes[0],data,constrictions,number_clusters,lambda);
    
    int first_wannabe;
    int second_wannabe;
    float prob = 0.7;

    //selección de los primeros contendientes
    do{
        first_wannabe = Randint(0, chromosomes.size() - 1);
        second_wannabe = Randint(0, chromosomes.size() - 1);
    } while(first_wannabe == second_wannabe);

    float first_wannabe_value = heuristicFunction(chromosomes[first_wannabe],data,constrictions,number_clusters,lambda);
    float second_wannabe_value = heuristicFunction(chromosomes[second_wannabe],data,constrictions,number_clusters,lambda);


    if(first_wannabe_value < second_wannabe_value)
    {
        selected_chromosomes[0] = chromosomes[first_wannabe];
    } else {
        selected_chromosomes[0] = chromosomes[second_wannabe];
    }

    // En el mismo bucle calculamos el mejor valor de los chromosomas precursores y hacemos los torneos binarios
    for(int i=1; i<chromosomes.size(); i++)
    {
        //torneo binario
        do
        {
            first_wannabe = Randint(0, chromosomes.size() - 1);
            second_wannabe = Randint(0, chromosomes.size() - 1);
        } while (first_wannabe == second_wannabe);

        first_wannabe_value = heuristicFunction(chromosomes[first_wannabe],data,constrictions,number_clusters,lambda);
        second_wannabe_value = heuristicFunction(chromosomes[second_wannabe],data,constrictions,number_clusters,lambda);

        if(first_wannabe_value < second_wannabe_value)
        {
            selected_chromosomes[i] = chromosomes[first_wannabe];
        }
        else
        {
            selected_chromosomes[i] = chromosomes[second_wannabe];
        }

        //busqueda del mejor valor
        float heuristic_value = heuristicFunction(chromosomes[i],data,constrictions,number_clusters,lambda);
        if(heuristic_value < best_chromosome_value){
            best_chromosome = chromosomes[i];
            best_chromosome_value = heuristic_value;
        }
    }

    FixedSegmentCrossoverOperator(selected_chromosomes, number_clusters, prob);
    //selected_chromosomes[worst_chromosome_index] = best_chromosome;

    //Aplicamos mutacion
    mutationOperator(selected_chromosomes, number_clusters);

    //Operador de arreglo
    for(int i = 0; i<selected_chromosomes.size(); i++){
        int isCorrect = correctSolution(selected_chromosomes[i],number_clusters);
        while(isCorrect != -1){
            selected_chromosomes[i][Randint(0,selected_chromosomes[i].size()-1)] = isCorrect;
            isCorrect = correctSolution(selected_chromosomes[i],number_clusters);
        }
    }


    int worst_position = 0;
    float worst_chromosome_value = heuristicFunction(selected_chromosomes[0],data,constrictions,number_clusters,lambda);
 
    for(int i = 1; i<selected_chromosomes.size(); i++){
        float heuristic_value = heuristicFunction(selected_chromosomes[i],data,constrictions,number_clusters,lambda);
        if(heuristic_value > worst_chromosome_value){
            worst_position = i;
            worst_chromosome_value = heuristic_value;
        }
    }

    selected_chromosomes[worst_position] = best_chromosome;

    heuristicEvaluations += chromosomes.size()*4; 

    return selected_chromosomes;
}

vector<vector<int>> UniformSegmentSteadySteateEvolutionScheme(const vector<vector<int>> &chromosomes, const vector<vector<float>> &data, const vector<vector<int>> &constrictions, int number_clusters, float lambda,  int &heuristicEvaluations)
{
    vector<vector<int>> chromosomes_copy = chromosomes;
    vector<vector<int>> selected_chromosomes(2);
    float prob = 1;

    for(int i=0; i<2; i++){
        int first_wannabe;
        int second_wannabe;

        //selección de los 2 padres
        do{
            first_wannabe = Randint(0, chromosomes.size() - 1);
            second_wannabe = Randint(0, chromosomes.size() - 1);
        } while(first_wannabe == second_wannabe);

        float first_wannabe_value = heuristicFunction(chromosomes[first_wannabe],data,constrictions,number_clusters,lambda);
        float second_wannabe_value = heuristicFunction(chromosomes[second_wannabe],data,constrictions,number_clusters,lambda);

        if(first_wannabe_value < second_wannabe_value)
        {
            selected_chromosomes[i] = chromosomes[first_wannabe];
        } else {
            selected_chromosomes[i] = chromosomes[second_wannabe];
        }
    }

    uniformCrossoverOperator(selected_chromosomes, number_clusters, prob);

    //Aplicamos mutacion
    mutationOperator(selected_chromosomes, number_clusters);

    //Operador de arreglo
    for(int i = 0; i<selected_chromosomes.size(); i++){
        int isCorrect = correctSolution(selected_chromosomes[i],number_clusters);
        while(isCorrect != -1){
            selected_chromosomes[i][Randint(0,selected_chromosomes[i].size()-1)] = isCorrect;
            isCorrect = correctSolution(selected_chromosomes[i],number_clusters);
        }
    }

    vector<float> values(chromosomes.size());
    for(int i = 0; i<chromosomes.size(); i++){
        values[i] = heuristicFunction(chromosomes[i],data,constrictions,number_clusters,lambda);
    }

    int worst_position = 0;
    float worst_chromosome_value = values[0];

    for(int i=1; i<values.size(); i++){
        if(values[i] > worst_chromosome_value){
            worst_position = i;
            worst_chromosome_value = values[i];
        }
    }

    values[worst_position] = 0;
    chromosomes_copy[worst_position] = selected_chromosomes[0];

    worst_position = 0;
    worst_chromosome_value = values[0];

    for(int i=1; i<values.size(); i++){
        if(values[i] > worst_chromosome_value){
            worst_position = i;
            worst_chromosome_value = values[i];
        }
    }

    values[worst_position] = 0;
    chromosomes_copy[worst_position] = selected_chromosomes[1];

    heuristicEvaluations += values.size()+4;

    return chromosomes_copy;
}

vector<vector<int>> FixedSegmentSteadySteateEvolutionScheme(const vector<vector<int>> &chromosomes, const vector<vector<float>> &data, const vector<vector<int>> &constrictions, int number_clusters, float lambda,  int &heuristicEvaluations)
{
    vector<vector<int>> chromosomes_copy = chromosomes;
    vector<vector<int>> selected_chromosomes(2);
    float prob = 1;

    for(int i=0; i<2; i++){
        int first_wannabe;
        int second_wannabe;

        //selección de los 2 padres
        do{
            first_wannabe = Randint(0, chromosomes.size() - 1);
            second_wannabe = Randint(0, chromosomes.size() - 1);
        } while(first_wannabe == second_wannabe);

        float first_wannabe_value = heuristicFunction(chromosomes[first_wannabe],data,constrictions,number_clusters,lambda);
        float second_wannabe_value = heuristicFunction(chromosomes[second_wannabe],data,constrictions,number_clusters,lambda);

        if(first_wannabe_value < second_wannabe_value)
        {
            selected_chromosomes[i] = chromosomes[first_wannabe];
        } else {
            selected_chromosomes[i] = chromosomes[second_wannabe];
        }
    }

    FixedSegmentCrossoverOperator(selected_chromosomes, number_clusters, prob);

    //Aplicamos mutacion
    mutationOperator(selected_chromosomes, number_clusters);

    //Operador de arreglo
    for(int i = 0; i<selected_chromosomes.size(); i++){
        int isCorrect = correctSolution(selected_chromosomes[i],number_clusters);
        while(isCorrect != -1){
            selected_chromosomes[i][Randint(0,selected_chromosomes[i].size()-1)] = isCorrect;
            isCorrect = correctSolution(selected_chromosomes[i],number_clusters);
        }
    }

    vector<float> values(chromosomes.size());
    for(int i = 0; i<chromosomes.size(); i++){
        values[i] = heuristicFunction(chromosomes[i],data,constrictions,number_clusters,lambda);
    }

    int worst_position = 0;
    float worst_chromosome_value = values[0];

    for(int i=1; i<values.size(); i++){
        if(values[i] > worst_chromosome_value){
            worst_position = i;
            worst_chromosome_value = values[i];
        }
    }

    values[worst_position] = 0;
    chromosomes_copy[worst_position] = selected_chromosomes[0];

    worst_position = 0;
    worst_chromosome_value = values[0];

    for(int i=1; i<values.size(); i++){
        if(values[i] > worst_chromosome_value){
            worst_position = i;
            worst_chromosome_value = values[i];
        }
    }

    values[worst_position] = 0;
    chromosomes_copy[worst_position] = selected_chromosomes[1];

    heuristicEvaluations += values.size()+4;

    return chromosomes_copy;
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
   // vector<int> solution;
    int number_clusters = stoi(argv[3]);
    int number_cromosomes = 50;
    int count;
    int seed = stoi(argv[4]);
    int change_count = 0;
    float lambda;

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

    cout << "\tAGG Segmento uniforme" << endl;
    
    auto start = chrono::steady_clock::now();
    int heuristicEvaluations = 0;
    count = 0;
    
    /*
     *  Generación de cromosomas
     */
    vector<vector<int>> chromosomes = generateIntitialChromosomes(data, number_clusters, number_cromosomes);

    /*
     * Aplicamos el esquema de evolución y de seleccion
     * 
     * */

    do{
        //chromosomes = FixedSegmentElitistEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
        chromosomes = UniformElitistEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
        //chromosomes = UniformSegmentSteadySteateEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);

    } while(heuristicEvaluations < 100000);

    /*
     * Escogemos la mejor solución 
     **/

    vector<int> best_chromosome = chromosomes[0];
    float best_chromosome_value = heuristicFunction(chromosomes[0],data,constrictions,number_clusters,lambda);

    for(int i = 1; i < chromosomes.size(); i++){
        float heuristic_value = heuristicFunction(chromosomes[i],data,constrictions,number_clusters,lambda);
        if(heuristic_value < best_chromosome_value){
            best_chromosome = chromosomes[i];
            best_chromosome_value = heuristic_value;
        }
    }
    auto end = chrono::steady_clock::now();

    // printSolution(best_chromosome, data);
    cout << "Tasa_C " << generalPartitionDesviation(best_chromosome, data, number_clusters) << endl;
    cout << "Tasa_inf " << infeasability(best_chromosome, constrictions) << endl;
    cout << "Agregado " << heuristicFunction(best_chromosome, data, constrictions, number_clusters, lambda) << endl;
    cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

    cout << "\tAGG Segmento fijo" << endl;
    
    start = chrono::steady_clock::now();
    heuristicEvaluations = 0;
    /*
     *  Generación de cromosomas
     */
    chromosomes = generateIntitialChromosomes(data, number_clusters, number_cromosomes);

    /*
     * Aplicamos el esquema de evolución y de seleccion
     * 
     * */

    do{
        chromosomes = FixedSegmentElitistEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
    } while(heuristicEvaluations < 100000);

    /*
     * Escogemos la mejor solución 
     **/

    best_chromosome = chromosomes[0];
    best_chromosome_value = heuristicFunction(chromosomes[0],data,constrictions,number_clusters,lambda);

    for(int i = 1; i < chromosomes.size(); i++){
        float heuristic_value = heuristicFunction(chromosomes[i],data,constrictions,number_clusters,lambda);
        if(heuristic_value < best_chromosome_value){
            best_chromosome = chromosomes[i];
            best_chromosome_value = heuristic_value;
        }
    }
    end = chrono::steady_clock::now();

    // printSolution(best_chromosome, data);
    cout << "Tasa_C " << generalPartitionDesviation(best_chromosome, data, number_clusters) << endl;
    cout << "Tasa_inf " << infeasability(best_chromosome, constrictions) << endl;
    cout << "Agregado " << heuristicFunction(best_chromosome, data, constrictions, number_clusters, lambda) << endl;
    cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

    cout << "\tAGE Segmento uniforme" << endl;
    start = chrono::steady_clock::now();
    heuristicEvaluations = 0;
    /*
     *  Generación de cromosomas
     */
    chromosomes = generateIntitialChromosomes(data, number_clusters, number_cromosomes);

    /*
     * Aplicamos el esquema de evolución y de seleccion
     * 
     * */

    do{
        chromosomes = UniformSegmentSteadySteateEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
    } while(heuristicEvaluations < 100000);

    /*
     * Escogemos la mejor solución 
     **/

    best_chromosome = chromosomes[0];
    best_chromosome_value = heuristicFunction(chromosomes[0],data,constrictions,number_clusters,lambda);

    for(int i = 1; i < chromosomes.size(); i++){
        float heuristic_value = heuristicFunction(chromosomes[i],data,constrictions,number_clusters,lambda);
        if(heuristic_value < best_chromosome_value){
            best_chromosome = chromosomes[i];
            best_chromosome_value = heuristic_value;
        }
    }
    end = chrono::steady_clock::now();

    // printSolution(best_chromosome, data);
    cout << "Tasa_C " << generalPartitionDesviation(best_chromosome, data, number_clusters) << endl;
    cout << "Tasa_inf " << infeasability(best_chromosome, constrictions) << endl;
    cout << "Agregado " << heuristicFunction(best_chromosome, data, constrictions, number_clusters, lambda) << endl;
    cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;

    cout << "\t AGE Segmento fijo" << endl;
    start = chrono::steady_clock::now();
    heuristicEvaluations = 0;
    /*
     *  Generación de cromosomas
     */
    chromosomes = generateIntitialChromosomes(data, number_clusters, number_cromosomes);

    /*
     * Aplicamos el esquema de evolución y de seleccion
     * 
     * */

    do{
        chromosomes = FixedSegmentSteadySteateEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
    } while(heuristicEvaluations < 100000);

    /*
     * Escogemos la mejor solución 
     **/

    best_chromosome = chromosomes[0];
    best_chromosome_value = heuristicFunction(chromosomes[0],data,constrictions,number_clusters,lambda);

    for(int i = 1; i < chromosomes.size(); i++){
        float heuristic_value = heuristicFunction(chromosomes[i],data,constrictions,number_clusters,lambda);
        if(heuristic_value < best_chromosome_value){
            best_chromosome = chromosomes[i];
            best_chromosome_value = heuristic_value;
        }
    }
    end = chrono::steady_clock::now();
    // printSolution(best_chromosome, data);
    cout << "Tasa_C " << generalPartitionDesviation(best_chromosome, data, number_clusters) << endl;
    cout << "Tasa_inf " << infeasability(best_chromosome, constrictions) << endl;
    cout << "Agregado " << heuristicFunction(best_chromosome, data, constrictions, number_clusters, lambda) << endl;
    cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;
    
    cout << "\tAM-(10,1.0) Fixed Segment Steady State + Local Search" << endl;
    
    start = chrono::steady_clock::now();
    heuristicEvaluations = 0;
    count = 0;
    
    /*
     *  Generación de cromosomas
     */
    chromosomes = generateIntitialChromosomes(data, number_clusters, number_cromosomes);

    int xi = chromosomes[0].size() * 0.1;
    int failures;
    vector<int> index;
    for (int i = 0; i < chromosomes[0].size(); i++)
    {
        index.push_back(i);
    }

    /*
     * Aplicamos el esquema de evolución y de seleccion
     * 
     * */

    do{
        //chromosomes = FixedSegmentElitistEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
        chromosomes = FixedSegmentSteadySteateEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
        //chromosomes = UniformSegmentSteadySteateEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
        count++;

        if(count == 10){
            random_shuffle(index.begin(), index.end());
            for(int sol_index=0; sol_index < chromosomes.size(); sol_index++){
                failures = 0;
                for (int i = 0; i < index.size() && failures < xi; i++)
                {
                    if (!lastElement(chromosomes[sol_index], index[i]))
                    {
                        int initial_cluster = chromosomes[sol_index][index[i]];
                        vector<int> partial_solution = chromosomes[sol_index];
                        partial_solution[index[i]] = 0;
                        float min_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                        int best_cluster = 0;
                        for (int j = 1; j < number_clusters; j++)
                        {
                            vector<int> partial_solution = chromosomes[sol_index];
                            partial_solution[index[i]] = j;
                            float function_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                            if (min_value > function_value)
                            {
                                min_value = function_value;
                                best_cluster = j;
                            }
                        }
                        if(best_cluster == initial_cluster){
                            failures++;
                        } else chromosomes[sol_index][index[i]] = best_cluster;
                        heuristicEvaluations += number_clusters;
                    }
                count = 0;
                }
            }
        }
        

    } while(heuristicEvaluations < 100000);

    /*
     * Escogemos la mejor solución 
     **/

    best_chromosome = chromosomes[0];
    best_chromosome_value = heuristicFunction(chromosomes[0],data,constrictions,number_clusters,lambda);

    for(int i = 1; i < chromosomes.size(); i++){
        float heuristic_value = heuristicFunction(chromosomes[i],data,constrictions,number_clusters,lambda);
        if(heuristic_value < best_chromosome_value){
            best_chromosome = chromosomes[i];
            best_chromosome_value = heuristic_value;
        }
    }
    end = chrono::steady_clock::now();

    // printSolution(best_chromosome, data);
    cout << "Tasa_C " << generalPartitionDesviation(best_chromosome, data, number_clusters) << endl;
    cout << "Tasa_inf " << infeasability(best_chromosome, constrictions) << endl;
    cout << "Agregado " << heuristicFunction(best_chromosome, data, constrictions, number_clusters, lambda) << endl;
    cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;
    index.clear();

    cout << "\tAM-(10,0.1) Fixed Segment Steady State + Local Search" << endl;
    
    start = chrono::steady_clock::now();
    heuristicEvaluations = 0;
    count = 0;
    
    /*
     *  Generación de cromosomas
     */
    chromosomes = generateIntitialChromosomes(data, number_clusters, number_cromosomes);

    xi = chromosomes[0].size() * 0.1;
    failures;
    //vector<int> index;
    for (int i = 0; i < chromosomes[0].size(); i++)
    {
        index.push_back(i);
    }

    /*
     * Aplicamos el esquema de evolución y de seleccion
     * 
     * */

    do{
        //chromosomes = FixedSegmentElitistEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
        chromosomes = FixedSegmentSteadySteateEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
        //chromosomes = UniformSegmentSteadySteateEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
        count++;

        if(count == 10){
            random_shuffle(index.begin(), index.end());
            for(int sol_index=0; sol_index < chromosomes.size(); sol_index++){
                if(Randfloat(0,1) < 0.1){
                    failures = 0;
                    for (int i = 0; i < index.size() && failures < xi; i++)
                    {
                        if (!lastElement(chromosomes[sol_index], index[i]))
                        {
                            int initial_cluster = chromosomes[sol_index][index[i]];
                            vector<int> partial_solution = chromosomes[sol_index];
                            partial_solution[index[i]] = 0;
                            float min_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                            int best_cluster = 0;
                            for (int j = 1; j < number_clusters; j++)
                            {
                                vector<int> partial_solution = chromosomes[sol_index];
                                partial_solution[index[i]] = j;
                                float function_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                                if (min_value > function_value)
                                {
                                    min_value = function_value;
                                    best_cluster = j;
                                }
                            }
                            if(best_cluster == initial_cluster){
                                failures++;
                            } else chromosomes[sol_index][index[i]] = best_cluster;
                            heuristicEvaluations += number_clusters;
                        }
                        count = 0;
                    }
                }
            }
        }
        

    } while(heuristicEvaluations < 100000);

    /*
     * Escogemos la mejor solución 
     **/

    best_chromosome = chromosomes[0];
    best_chromosome_value = heuristicFunction(chromosomes[0],data,constrictions,number_clusters,lambda);

    for(int i = 1; i < chromosomes.size(); i++){
        float heuristic_value = heuristicFunction(chromosomes[i],data,constrictions,number_clusters,lambda);
        if(heuristic_value < best_chromosome_value){
            best_chromosome = chromosomes[i];
            best_chromosome_value = heuristic_value;
        }
    }
    end = chrono::steady_clock::now();

    // printSolution(best_chromosome, data);
    cout << "Tasa_C " << generalPartitionDesviation(best_chromosome, data, number_clusters) << endl;
    cout << "Tasa_inf " << infeasability(best_chromosome, constrictions) << endl;
    cout << "Agregado " << heuristicFunction(best_chromosome, data, constrictions, number_clusters, lambda) << endl;
    cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;
    index.clear();

    cout << "\tAM-(10,0.1mej) Fixed Segment Steady State + Local Search" << endl;
    
    start = chrono::steady_clock::now();
    heuristicEvaluations = 0;
    count = 0;
    /*
     *  Generación de cromosomas
     */
    chromosomes = generateIntitialChromosomes(data, number_clusters, number_cromosomes);

    xi = chromosomes[0].size() * 0.1;
    int best_chromosomes = chromosomes.size() * 0.1;
    failures;
    //vector<int> index;
    for (int i = 0; i < chromosomes[0].size(); i++)
    {
        index.push_back(i);
    }

    /*
     * Aplicamos el esquema de evolución y de seleccion
     * 
     * */

    do{
        //chromosomes = FixedSegmentElitistEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
        chromosomes = FixedSegmentSteadySteateEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
        //chromosomes = UniformSegmentSteadySteateEvolutionScheme(chromosomes,data,constrictions,number_clusters,lambda, heuristicEvaluations);
        count++;

        if(count == 10){
            vector<pair<float, int>> sorted_values;
            for(int i = 0; i < chromosomes.size(); i++){
                sorted_values.push_back( make_pair(heuristicFunction(chromosomes[i], data, constrictions, number_clusters, lambda),i) );
            }
            heuristicEvaluations += chromosomes.size();

            sort(sorted_values.begin(), sorted_values.end());

            for(int sol_index=0; sol_index < best_chromosomes; sol_index++){
                random_shuffle(index.begin(), index.end());
                    failures = 0;
                    for (int i = 0; i < index.size() && failures < xi; i++)
                    {
                        if (!lastElement(chromosomes[sorted_values[sol_index].second], index[i]))
                        {
                            int initial_cluster = chromosomes[sorted_values[sol_index].second][index[i]];
                            vector<int> partial_solution = chromosomes[sorted_values[sol_index].second];
                            partial_solution[index[i]] = 0;
                            float min_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                            int best_cluster = 0;
                            for (int j = 1; j < number_clusters; j++)
                            {
                                vector<int> partial_solution = chromosomes[sorted_values[sol_index].second];
                                partial_solution[index[i]] = j;
                                float function_value = heuristicFunction(partial_solution, data, constrictions, number_clusters, lambda);
                                if (min_value > function_value)
                                {
                                    min_value = function_value;
                                    best_cluster = j;
                                }
                            }
                            if(best_cluster == initial_cluster){
                                failures++;
                            } else chromosomes[sorted_values[sol_index].second][index[i]] = best_cluster;
                            heuristicEvaluations += number_clusters;
                        }
                        count = 0;
                    }
                
            }
        }
        

    } while(heuristicEvaluations < 100000);

    /*
     * Escogemos la mejor solución 
     **/

    best_chromosome = chromosomes[0];
    best_chromosome_value = heuristicFunction(chromosomes[0],data,constrictions,number_clusters,lambda);

    for(int i = 1; i < chromosomes.size(); i++){
        float heuristic_value = heuristicFunction(chromosomes[i],data,constrictions,number_clusters,lambda);
        if(heuristic_value < best_chromosome_value){
            best_chromosome = chromosomes[i];
            best_chromosome_value = heuristic_value;
        }
    }
    end = chrono::steady_clock::now();

    // printSolution(best_chromosome, data);
    cout << "Tasa_C " << generalPartitionDesviation(best_chromosome, data, number_clusters) << endl;
    cout << "Tasa_inf " << infeasability(best_chromosome, constrictions) << endl;
    cout << "Agregado " << heuristicFunction(best_chromosome, data, constrictions, number_clusters, lambda) << endl;
    cout << "Tiempo " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " ms" << endl;
    index.clear();
}