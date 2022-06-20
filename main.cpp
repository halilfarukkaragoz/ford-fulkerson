/**
 * @file main.cpp main file of the solution // Only have main function which calls final functions
 * @author Halil Faruk Karagöz 150180014  (karagozh18@itu.edu.tr)
 * @brief Solutions of third Assignment of the AA2 lecture
 * Problem : 
 * There are 2 input : Adjacency matrix of a graph and a matrix which inform possibilities of closing edges
 * Find the maximum flows with the possibilities of happening 
 * @date 2022-06-02
 * @copyright Copyright (c) 2022
 */

#include "graph.cpp"


int main(int argc,  char *argv[]){

    /* 
        Argument taken from terminal
    */

    string file_name = argv[1];
    int start_p = stoi(string(argv[2]));
    int end_p = stoi(string(argv[3]));
        
    /* 
        Read file and fill 2 matrices 
        matrix1 = Adjacent matrix of given graph 
        matrix2 = Possibilities of edges to be closed 
    */
    d2array matrix1;
    d2array matrix2;  
    read_input(file_name,matrix1,matrix2);
    
    /* Instance of graph */
    Graph graph(matrix1);

    /* 
        Generate possible graphs by the possibilty matrix and find maximum flow of them 
     */
    vector<pair<float,float> > v = generate_possiblity_flow_pairs(graph,matrix2,start_p,end_p);

    /* 
        Print output as desired  
    */
    print_output(v,start_p,end_p);

    return 0;
}
