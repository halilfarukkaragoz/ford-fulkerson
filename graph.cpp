/**
 * @file graph.cpp Have graph related Functions
 * @author Halil Faruk Karagöz 150180014  (karagozh18@itu.edu.tr)
 * @version 0.1
 * @date 2022-06-02
 * @copyright Copyright (c) 2022
 */
#include "graph.hpp"
#include <string>
#include <fstream>
#include<vector>
#include<algorithm>
#include <math.h>
#include <map>

/**
 * @brief Create a 2dmatrix object
 * 
 * @param row number of row
 * @param col number of column
 * @param matrix matrix which is to be resized
 */
void create_2dmatrix(int row,int col,d2array &matrix){
    matrix.resize(row);

    for(int i = 0; i < row; i++){
        matrix[i].resize(col);
    }
}
void print_matrix(d2array adjacent_matrix){
        for(int i =0; i < adjacent_matrix.size();i++){ 
            for(int j = 0; j < adjacent_matrix[0].size(); j++){
                cout<<adjacent_matrix[i][j]<< " ";
                }
            cout<<endl;  
        }  
}

/**
 * @brief It read inout from file and return a 2d array
 * 
 * @param file_name name of the file which input will be read
 */
void read_input(string file_name,d2array &matrix1,d2array &matrix2){
    ifstream read_file(file_name);

    int node_number = 0;
    string row,word;
    if(read_file.is_open()){

        getline(read_file,row);
        node_number = stoi(row);
        create_2dmatrix(node_number,node_number,matrix1);
        create_2dmatrix(node_number,node_number,matrix2);
        for(int i = 0; i < node_number; i++){
            for(int j = 0; j < node_number; j++){
                getline(read_file,word,' ');
                matrix1[i][j] = stoi(word);
            }
        }

    }


    for(int i = 0; i < node_number; i++){
        for(int j = 0; j < node_number; j++){
            getline(read_file,row,' ');
            matrix2[i][j] = stof(row);
        }
    }
}

vector<int> trace_back(vector<int> &route,int starting_p,int ending_p){
    vector<int> path;
    int current_node = ending_p-1;
    path.push_back(current_node);
    while(current_node != starting_p-1){
        current_node = route[current_node];
        path.push_back(current_node);
    }
    reverse(path.begin(),path.end());
    return path;
}

/**
 * @brief augment a path using BFS
 *
 * @param g graph to be traversed 
 * @return vector<int>: vector which for augmented path 
 * -----------------
 * @complexity: 
 * In worst case scenario: Traverse all edges and all nodes in a graph O(N+E)
 */
vector<int> augment_path(Graph &g,int start_p,int end_p){
    /* Queue for adding next nodes to traverse */
    vector<int> q; 
    /* Keep traversed nodes to prevent go same node twice */
    vector <int> traversed_nodes; 
    /* Adjacency matrix of graph */
    d2array adjacent_matrix = g.get_adjacency_mat(); 
    /* Keep the where is it come from to trace back route */
    vector <int> route(adjacent_matrix.size(),0);
    /* Push starting node as init node */
    q.push_back(start_p-1);
    /* As long as queue is empty keep continue to search*/
    while(q.size() != 0){
        /* Take next node to traverse from queue, erase it from the queue and add it to the vector that holds traversed nodes*/
        int current_node = q.front();
        q.erase(q.begin());
        traversed_nodes.push_back(current_node);
        /* Check edges in current node */
        for(int i = 0; i <adjacent_matrix.size(); i++){
            if(adjacent_matrix[current_node][i] != 0 ){
                /* If there is an edge that current node to ending node */ 
                if(i == end_p-1){
                    /* Push ending node to vector which hold traversed nodes */ 
                    traversed_nodes.push_back(i);
                    /* Put the information of where did it come from to the last node to trace back route */
                    route[i] = current_node;
                    /* Trace back to route from ending point to starting point */
                    return trace_back(route,start_p,end_p);
                }
                else {
                    /* Check is node is already traversed or already in the queue */
                    bool add_queue = true;
                    for(int j = 0; j < traversed_nodes.size(); j++)
                        if(traversed_nodes[j] == i){
                            add_queue = false;
                            break;
                        }
                    for(int j = 0; j < q.size(); j++)
                        if(q[j] == i){
                            add_queue = false;
                            break;
                        }    
                    /* If it is not already traversed or already in the queue push element to the queue to traverse later */
                    if(add_queue){ 
                        q.push_back(i);
                        route[i] = current_node;
                    }
                }    
            }
        }
    }
    /* If there is no path return an empty vector */
    vector<int> empty_vector;
    return empty_vector;
}


/**
 * @brief Update residual graph by the augmented path 
 * 
 * @param g graph to be updated 
 * @param path path to be updated 
 * @return int: 
 * ---------------
 * @complexity: 
 * In worst case path include all edges so iterate edge number O(E)
 */
int update_graph(Graph &g,vector<int> path){
    /* Get adjacent matrix of current graph */
    d2array ad_matrix =g.get_adjacency_mat();
    /* To find bottleneck of the path */
    int min = ad_matrix[path[0]][path[1]];
    /* Check all edges in the path and compare to find bottleneck */
    for(int i = 0; i < path.size()-1; i++){
        if (ad_matrix[path[i]][path[i+1]] < min) min = ad_matrix[path[i]][path[i+1]];
    }
    /* Update all edges with the found min element */
    for(int i = 0; i < path.size()-1; i++){
        ad_matrix[path[i]][path[i+1]] -= min;
        ad_matrix[path[i+1]][path[i]] += min;
    }
    g.assign_adjacency_mat(ad_matrix);
    return min;
}


/**
 * @brief find max flow of given graph
 * 
 * @param g 
 * @return int: max flow value  
 * ----------------------
 * @complexity: 
 * In worst case augment path works in O(n+e) 
 * Possible paths * O(n+e) 
 */
int ford_fulkerson(Graph g,int start_p,int end_p){
    /* Residual Graph for updating depends on the augmented path */
    Graph res_g = g;
    /* Flag for stoping condition */
    bool stop = false;
    /* Keep the sum of flow */
    int max_flow = 0;
    /* Loop run unless there is no path between start point and end point  */
    while(!stop){
        /* Augment a path */
        vector<int> path = augment_path(res_g,start_p,end_p);
        /* If path size is 0 then stop the algorithm */
        if(path.size() == 0){
            stop = true;
            break;
        }
        /* Update residual graph and add flow of the current path to sum of the flow */
        max_flow += update_graph(res_g,path);
    }
    /* Return sum of the flows */
    return max_flow;
}

/**
 * @brief Return row and column information of nonzero element
 * 
 * @param G Graph to be reviewed 
 * @return d2array : return row and columns of all nonzero elements in 2d matrix
 */
d2array return_nonzero_elements(d2array mat){
    d2array nonzero_elements;
    for(int i = 0; i < mat.size();i++){
        for(int j  = 0; j < mat.size(); j++){
            if (mat[i][j] != 0){
                d1array v;
                v.push_back(i);v.push_back(j);nonzero_elements.push_back(v); // Push row, and column to a 1d vector than push 1d vector to 2dvector
            }
        }
    }
    return nonzero_elements;
}

/**
 * @brief Calculate possibility of current situation 
 * 
 * @param possibility_matrix Possibility matrix that given as input 
 * @param nonzer_edges Places of nonzero edges (x,y)
 * @param current Current closed edge 
 * @return float : possibility of current situation 
 */
float calculate_poss(d2array possibility_matrix,d2array nonzer_edges,vector<bool> current){
    float poss = 1;
    for(int i = 0; i < nonzer_edges.size();i++){
        int row = nonzer_edges[i][0];int col = nonzer_edges[i][1];
        if(current[i])
            poss *=possibility_matrix[row][col];
        else
            poss *=(1 - possibility_matrix[row][col]);

    }
    return poss;

}

/**
 * @brief To sort vector<pair<float,float> descending order
 * It will be used in std::sort() function 
 * 
 * @param a first element to compare 
 * @param b second element to compare 
 * @return true when first element > second 
 * @return false when first element < second 
 */
bool cmp(pair<float, float>& a,
         pair<float, float>& b)
{
    return a.second > b.second;
}

/**
 * @brief Generate all possible graphs, then find probability of that graph to be happen with maximum flow of that graph 
 * 
 * @param G  Graph to be calculated max_flow - probability duo 
 * @param possibilty_matrix Possibility of edges to be closed
 * @param start_p Starting node of flow
 * @param end_p Ending node of flow
 * @return vector<pair<float,float> > Keeps max_flow - probability pairs in descending order by probabilty  
 */
vector<pair<float,float> > generate_possiblity_flow_pairs(Graph G, d2array possibilty_matrix,int start_p,int end_p){
    /* Used map data structures to keep (key : FLow , value : Possiblity) structure */
    map<float,float> possibilities_dict;

    /* 
        Return position(x,y) of the edges which has nonzero possibilties 
        @complexity: Traverse through all element in adjacent matrix O(n^2) : n = node number in graph 
    */
    d2array nonzero_poss = return_nonzero_elements(possibilty_matrix);
    
    /* Number of subset */ 
    int it_number = pow(2,nonzero_poss.size());

    /* Sum of all possibilties is 1 at the beginning */
    float remaining_poss = 1;

    /* 
        Generate possible graph for all subset then find flow for generated graph 
        Keep results in map
        ---------------------
        @complexity: In worst set scenario, between all nodes there is an edge and all edges have possibility to closed 
        node number = n
        edge number : n * (n-1)
        subset : 2^(n^2) 
        For loop iterate 2^(n^2) times. O(2^(n^2)) :(
    */
    for(int i = 0; i < it_number; i++){
        /* To keep current closed edges */
        vector<bool> current_array(nonzero_poss.size(),false);
        /* Pick edges to close in order */
        for(int j = 0; j< nonzero_poss.size();j++){
            if (i & (1 << j)){
                current_array[j] = true;
            }
        }
        
        /* Calculate possibility of subset to happen */
        float poss = calculate_poss(possibilty_matrix,nonzero_poss,current_array);
        
        /* Generate adjacent matrix of current subset */
        d2array adjacent_mat = G.get_adjacency_mat();
        for(int j = 0; j <nonzero_poss.size();j++){
            int row = nonzero_poss[j][0];int col = nonzero_poss[j][1];
            if(current_array[j])
                adjacent_mat[row][col] = 0;
        }

        if (i == 0) continue; // empty set case this will be calculated at the end
        
        /* Generate graph for the put in Ford-Fulkerson Algorithm, then run the ford fulkerson algorithm*/
        Graph g_new(adjacent_mat);
        int flow = ford_fulkerson(g_new,start_p,end_p);
        /* If flow is already exist add increase possibility with the current possibility */
        if(possibilities_dict.find(flow) != possibilities_dict.end()){
            possibilities_dict[flow] +=poss;
        }
        /* If it doesn't exist put a new flow-poss pair to map data structures */
        else possibilities_dict.insert(pair<float,float>(flow,poss));
        remaining_poss -=poss;
    }

    /* Empty subset case, possibilty = remaining possibilty */
    int flow = ford_fulkerson(G,start_p,end_p);
    if(possibilities_dict.find(flow) != possibilities_dict.end()){
            possibilities_dict[flow] +=remaining_poss;
        }
    else possibilities_dict.insert(pair<float,float>(flow,remaining_poss));
    
    /* 
        To sort possibilities in map, at first put element to a vector 
        Then used std::sort function to sort
    */
    vector<pair<float,float> > possiblity_flow_pairs;
        for(std::map<float,float>::iterator it = possibilities_dict.begin(); it != possibilities_dict.end(); ++it) {
        possiblity_flow_pairs.push_back(pair<float,float>(it->first,it->second));
    }
    sort(possiblity_flow_pairs.begin(),possiblity_flow_pairs.end(),cmp);

    return possiblity_flow_pairs;
}


/**
 * @brief 
 * 
 * @param v 
 * @param start_p 
 * @param end_p 
 * @result Printed flow-probability duo as desired
 */
void print_output(vector<pair<float,float> > v, int start_p, int end_p){
    /*  
        Calculate expected flow 
        calculated flow = flow_1 * poss_1 + flow_2 * poss_2 .... flow_n * poss_n (In all element)
    */
    float expected_flow = 0;

    /* Print output in desired way */
    for(vector<pair<float,float> >::iterator it = v.begin(); it != v.end(); ++it) {
    expected_flow +=(it->first * it->second);
    std::cout << "Probability of occurence:  " << it->second << ", ";
    std::cout << "Maximum Flow: " << it->first <<endl;
    }

    cout <<"Expected Maximum Flow from node "<<start_p <<" to node "<<end_p<<": "<<expected_flow<<endl;


}