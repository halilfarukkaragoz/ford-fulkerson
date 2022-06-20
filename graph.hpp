/**
 * @file graph.hpp Definition of Graph class 
 * @author Halil Faruk Karagöz 150180014  (karagozh18@itu.edu.tr)
 * @version 0.1
 * @date 2022-06-02
 * @copyright Copyright (c) 2022
 */

#pragma once 

#include <iostream>
#include <vector>
using namespace std;
#define d1array vector<float>
#define d2array vector<d1array>


/**
 * @brief A graph representation with adjacent matrix 
 */
class Graph{
    d2array adjacent_matrix; // Edge values show the capacity 
    public:
    Graph(d2array matrix1){adjacent_matrix = matrix1;}
    void assign_adjacency_mat(d2array matrix1){adjacent_matrix = matrix1;}
    d2array get_adjacency_mat(){return adjacent_matrix;}
};




