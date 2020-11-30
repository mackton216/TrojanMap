# EE599 Final Project Report
## Travelling Trojan Problem
### Ashwin Sivakumar, Mackton Vishal Mendonca

## Overview of Design


<p align="left"><img src="img/TravellingSalesman.png" alt="TSP" width="500"/></p>

- The travelling salesman problem is a well known [NP-hard problem](https://en.wikipedia.org/wiki/NP-hardness) that is used in computer science for combinatorial optimizations.
- Our design employs the use of structures and helper functions to try and solve this problem.
- We have made the use of the Haversine formula to find the distance in miles between various points differing in longitude and latitudes.
- We use the Dijkstra algorithm along with a priority queue to find the Shortest Path between 2 points on the map.
- We use Brute Force and several helper functions to reach our goal in helping our Travelling Trojan reach all the nodes in a given map and return in the shortest path possible.
- We have also used two heuristic methods to find “good enough outputs” for the Travelling Salesman Problem namely:
 - 2 - Opt
 - ##### Genetic Algorithm (extra credit implementation)
## The Node Class
```cpp
class Node {
 public:
 Node(){};
 Node(const Node &n){id = n.id; lat = n.lat; lon = n.lon; name = n.name; neighbors = n.neighbors;};
 std::string id; // A unique id assign to each point
 double lat; // Latitude
 double lon; // Longitude
 std::string name; // Name of the location. E.g. "Bank of America".
 std::vector<std::string>neighbors; // List of the ids of all neighbor points.
};
```
- Contains a parameterized constructor with datatype Node
- The constructor defines objects for Node IDs, latitudes, longitudes, names of places represented by the nodes and neighbouring nodes

## Primary Code Functions:
```cpp
// Returns a vector of names given a partial name.
std::vector<std::string> Autocomplete(std::string name);
// Given the name of two locations, it should return the **ids** of the nodes on the shortest path using Dijkstra Priority Queue
std::vector<std::string> CalculateShortestPath(std::string location1_name, std::string location2_name);
// Travelling Salesman Using Brute Force
std::pair<double, std::vector<std::vector<std::string>>> TravellingTrojan(std::vector<std::string> &location_ids);
// Travelling Salesman Using 2-Opt
std::pair<double, std::vector<std::vector<std::string>>> TravellingTrojan_2opt(std::vector<std::string> &location_ids);
// Travelling Salesman Using Genetic Algorithm
std::pair<double, std::vector<std::vector<std::string>>> TravellingTrojan_GeneticAlgo(std::vector<std::string> &location_ids);
```
- **Autocomplete** - To take a part of a location name from the user and return all places that contain that part of the name.
- **GetPosition** - To find the latitude and longitude of a given location
- **CalculateShortestPath** - Calculate the shortest path between 2 locations
- **TravellingTrojan** - Given a set of locations, find the shortest route between after visiting every location once and returning using Brute Force.
- **Heuristical solutions to TravellingTrojan using 2_opt and GeneticAlgo.**

## Updated Menu
<p align="center"><img src="img/NewMenu.PNG" alt="NewMenu" width="500"/></p>

- Our updated menu allowing 7 options.

## Generating the weight matrix
<p align="center"><img src="img/GenWeightMatrix.PNG" alt="GenWeightMatrix" width="500"/></p>

- Objective was to create an adjacency matrix while mapping every location id to an index.
- This way the core functions needed minimal alteration.
- Made implementation more modular.

## Shortest Path - Dijkstra Priority Queue
<p align="center"><img src="img/ShortestPathDjikstra.PNG" alt="SP_Dijkstra" width="500"/></p>

- Idea is to traverse through graph and find next minimum distanced node and update distances of all children.
- Store parent nodes for each node in every iteration.
- Leverages Priority queue functionality to efficiently obtain next minimum unvisited node.
- Parent nodes are accessed and recursed to obtain the path in addition to the shortest distance.

## Travelling Salesman – Brute Force
<p align="center"><img src="img/TSPBruteForce.PNG" alt="TSP_brute" width="500"/></p>

- Aimed at solving the TSP using a greedy technique. 
- Method always returns an optimal result but computationally expensive.
- Intuitive technique inspired from DFS search + Permutations.

## Travelling Salesman – 2 opt
<p align="center"><img src="img/TSP2-Opt.PNG" alt="TSP_2Opt" width="500"/></p>

- Aimed at solving the TSP using an approximate technique.
- Utilizes 2-edge swapping
- Method does not always return optimal distance / path , but it is expected to give a “very good” approximate of the ideal output.
- Computationally less expensive in comparison with Brute Force.

## Travelling Salesman – Genetic Algorithm
<p align="center"><img src="img/GA_FC.png" alt="TSP_GA" width="500"/></p>


- Initialize the population randomly.
- Determine the fitness (path length) of the chromosome (string containing path progression).
- Until done repeat:
 - Select parents.
 - Perform crossover and mutation.
 - Calculate the fitness of the new population.
 - Append it to the gene pool.

---

## Every function explained:
### GetLat, GetLon, GetName:
```cpp
double TrojanMap::GetLat(std::string id) 
{ 
 return data[id].lat;
}

double TrojanMap::GetLon(std::string id) 
{ 
 return data[id].lon;
}

std::string TrojanMap::GetName(std::string id) 
{ 
 return data[id].name;
}
```
- Using the attributes of the variable ‘data’ defined in the cpp.h file.
- Data is a map consisting of strings and Nodes defined in the .csv file.
- String “id” is user input and is mapped to the latitude, longitude and name (if any) of the respective Node id.
- Each of these have a time complexity of O(1).

### GetNeighborIDs
```cpp
std::vector<std::string> TrojanMap::GetNeighborIDs(std::string id) 
{ 

 std::vector<std::string> result;

 for(int i = 0; i < data[id].neighbors.size(); i++)
 {
 result.push_back(data[id].neighbors[i]);
 }
 return result;
} 
```
- Uses the neighbors attribute of the map ‘data’ in the cpp.h file to retrieve the neighboring IDs of each Node ID as defined the .CSV file.
- Returns a vector of strings.
- O(N) where N is the number of nodes

### CalculateDistance
```cpp
 // TODO: Use Haversine Formula:
 // dlon = lon2 - lon1;
 // dlat = lat2 - lat1;
 // a = (sin(dlat / 2)) ^ 2 + cos(lat1) * cos(lat2) * (sin(dlon / 2)) ^ 2;
 // c = 2 * arcsin(min(1, sqrt(a)));
 // distances = 3961 * c;

 // where 3961 is the approximate radius of the earth at the latitude of
 // Washington, D.C., in miles
double TrojanMap::CalculateDistance(const Node &a, const Node &b) 
{

 double dlon, dlat, a1, c, distance;

 dlon = b.lon - a.lon;
 dlon = dlon * (M_PI/180);
 dlat = b.lat - a.lat;
 dlat = dlat * (M_PI/180);

 a1 = (pow((std::sin(dlat / 2)),2)) + std::cos(a.lat* (M_PI/180)) * std::cos(b.lat* (M_PI/180)) * (pow((std::sin(dlon / 2)),2));
 c = 2 * std::asin(MIN(1, sqrt(a1)));

 distance = (3961 * c);

 return distance;
}
```
- Using the Haversine formula to calculate distance between 2 Nodes which contain respective latitudes and longitudes. 
- Returns a double.
- O(1)

### CalculatePathLength
```cpp
double TrojanMap::CalculatePathLength(const std::vector<std::string> &path) 
{ 

 int i = 0, j = 1;
 double sum = 0;
 while(j<path.size())
 {
 sum += CalculateDistance(data[path[i]], data[path[j]]);
 i++;
 j++; 
 }
 
 return sum;
}
```
- Uses CalculateDistance function to find distance between the latitudes and longitudes of neighboring nodes that are in the path.
- Sums it up to find total length of a path taken.
- Returns a double.
- O(n).

### GetNode and GetNodeFromName
```cpp
Node TrojanMap::GetNode(std::string id)
{
 std::map<std::string, Node>::iterator it;

 for(it = data.begin(); it!=data.end();it++)
 {
 if(it->first.compare(id) == 0)
 {
 return it->second;
 break; 
 }
 }
}


Node TrojanMap::GetNodeFromName(std::string locname)
{
 std::map<std::string, Node>::iterator it;
 int flag = 0;
 for(it = data.begin(); it!=data.end();it++)
 {
 if(it->second.name.compare(locname) == 0)
 {
 flag = 1;
 return it->second;
 break; 
 }
 }
}
```
- The id of the node is compared to the ids present in the .CSV file and the corresponding node is returned.
- O(n)
- The location name is searched for in the .CSV file and then the corresponding Node is returned. 
- O(n)

### Autocomplete
```cpp
std::vector<std::string> TrojanMap::Autocomplete(std::string name) 
```
- An input from the user is stored as a string and converted into lower case.
- A map is created with the string as the key and the Node as the data.
- A temporary string is created to hold a part of the name. This temp string is used to scan the .CSV file and return strings which match the temp string.
- If they match the temp string, the string is pushed into the results vector and returned to the user.
- O(n)

### GetPosition
```cpp
std::pair<double, double> TrojanMap::GetPosition(std::string name) {
 std::pair<double, double> results(-1, -1);
 std::map<std::string, Node>::iterator it;

 for(it = data.begin(); it!=data.end();it++)
 {
 if(it->second.name.compare(name) == 0)
 {
 results.first = it->second.lat;
 results.second = it->second.lon;
 return results;
 break; 
 }
 }
}
```
- An iterator is used to traverse the .CSV file holding the nodes.
- This iterator compares the name string that the user inputs with the name string in the file.
- If there is a match in the name, the corresponding latitudes and longitudes of the name are retrieved from the .CSV file, stored in a pair and displayed to the user.
- O(n)

## Calculate the Shortest Path - A single source shortest path approach.
<p align="center"><img src="img/Dijkstra_Animation.gif" alt="TSP_GA" width="300"/></p>

### DikjstraPriorityQueue

```cpp
std::vector<double> TrojanMap::DijkstraPriorityQueue(int source, std::vector<std::vector<double>> weight_, std::map<int, int> &prev)
```
- Uses the Dijkstra algorithm to find the shortest path between all the nodes in the map.
- The priority queue is used to make it more efficient to find the shortest path.
- The function starts at the source node.
- The priority queue is used to extract the next closest unvisited node from the list of currently visited ones.
- Once found the distances of all its neighbors are updated.
- This continues till all nodes in the graph are visited.
- Ultimately returns the shortest path taken to cover all the nodes.
- We use this function to find the shortest path taken between 2 nodes defined by the user. 
- **Time Complexity** O((m+n)logn)

### CalculateShortestPath

```cpp
std::vector<std::string> TrojanMap::CalculateShortestPath(std::string location1_name, std::string location2_name) 
```
- This function takes in 2 strings from the user. These strings represent the names of the places on the map.
- The GetNodeFromName function then retrieves the nodes with their respective latitudes and longitudes from the .CSV file and they are assigned to 2 integer values.
- Starting from the source node, the distances of the neighbouring nodes are calculated. 
- The neighbouring node that is the closest to the source is found and the distance between them is stored. This is repeated for all the nodes connecting the source to the destination node. 
- A weight matrix is created from evaluation of distances of neighbouring nodes by pushing back these distances into a vector of double vectors.
- The reverse of the path generated from the source to destination is evaluated to verify the path.
- These values are passed as parameters for the DijkstraPriorityQueue function.
- The path is returned.

- **Time Complexity** O(n^2).

## Travelling Salesman .. a) A Brute Force Approach

### TSPHELPER

```cpp
double TrojanMap::TSPHELPER(int start, int cur_node, double cur_cost, std::vector<std::string> path, std::vector<std::vector<double>> weight_, std::map<int, std::string>& rev_index, std::pair<double, std::vector<std::vector<std::string>>>&progress)
```
- The TSPHELPER function is a function written to evaluate the various paths that the Trojan can take to reach a destination and come back to the source node.
- The parameters passed here are the start int, the current node, the “cost” of travelling through the nodes, the path taken, the a vector of vectors of weights, a map of integers and strings containing the indices in reverse and a pair of double and vector of vector of strings containing the total progress.
- Each node that the Trojan visits, the distance between the nodes are evaluated and stored.
- Each path that the Trojan visits is stored in a vector.
- This happens recursively for each of the neighbouring nodes from the source to the destination.
- Once all the paths have been traversed, the path with the least distance is chosen.

### TravellingTrojan

```cpp
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan(std::vector<std::string> &location_ids) 
```
- The parameters passed are a vector of strings containing the location ids.
- The function is used to create a weight matrix. This weight matrix contains the various distances covered by the Trojan while traversing through the map
- These sizes are passed to the CalculateDistance function to find the total distance the Trojan travels after passing through each node between source and destination once. 
- The progress i.e. node ids in every path chosen are stored. 
- The weight matrix along with the first node, final node, paths taken, the reverse path and the progress of each path is passed to the TSPHELPER function to pick the shortest path.
- The CalculatePathLength function is called to find the total distance of the shortest path in miles.
- The various nodes covered in path determined as the shortest are returned.
- Returns the shortest path length and the progression vectors for the animation. 
- **Time Complexity** O(n!)


## Travelling Salesman .. b) 2 opt Heuristic Approach


<p align="center"><img src="img/2-OptImage.PNG" alt="algo_image" width="300"/></p>

<p align="center"><img src="img/2-opt.gif" alt="algo_image" width="300"/></p>

**Note: Interpretation of output while executed from the menu**

<p align="center"><img src="img/2_optOP.png" alt="2_opt_OP" width="500"/></p>

- Once the input number of locations is given, the ouput will consist of initial, final and difference in path lengths.
- This happens before the shortest path is displayed. (Easy to miss the summary mentioned in the first point)
- Tests cases pass if an improvement is made.

### NN
```cpp
std::vector<int> TrojanMap::NN(std::vector<std::vector<double>> weight_, double &dist)
```
- This is the function to obtain the nearest neighbor heuristic, as a starting point for the 2-opt algorithm
- Execution starts of at a node, and the algorithm identifies the next closest node and establishes a link.
- Now the current node is the newly found point and the previous step is repeated till all nodes are covered.

### Distance
```cpp
double TrojanMap::distance(std::vector<int> pt, int i, std::vector<std::vector<double>> weight_)
```
- Calculates the distance of a cyclic vector during 2-opt execution.

### TWO_OPT

```cpp
std::pair<double, std::vector<std::vector<int>>> TrojanMap::TWO_OPT(std::vector<std::vector<double>> weight_, std::pair<double, std::vector<std::vector<int>>> result)
```
- This is the core part of the 2-opt algorithm
- Takes in the the nearest neighbor heuristic as its initial input.
- Identifies possible edge swaps with the current path and implements all possible edge swap.
- If an improvement in the path lenght is found, the loop is broken and we start over swapping edges with the new shortest path.
- The algorithm continues till no improvement is made after all edge swaps are complete.

### TravellingTrojan_2opt
```cpp
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan_2opt(std::vector<std::string> &location_ids)
```
- The function generates the weight matrix based on input location ids.
- Invokes the nearest neighbor heuristic and then the 2opt algorithm.
- Returns the shortest path length and the progression vectors for the animation. 


## Travelling Salesman .. c) Genetic Algorithm

- Tests cases pass if an improvement is made.

<p align="center"><img src="img/GenAlgo.GIF" alt="GA_op" width="500"/></p>

**Note: Interpretation of output while executed from the menu**

<p align="center"><img src="img/GA_OP.png" alt="GA_op" width="500"/></p>

```cpp
int rand_num(int start, int end) 
```
- Function to return a random number from start and end.

```cpp
bool repeat(string s, char ch) 
```
- Function to check if the character has already occurred in the string.

```cpp
string mutatedGene(string gnome)
```
- Function to return a mutated GNOME Mutated GNOME is a string with a random interchange of two genes to create variation in species. 

```cpp
string create_gnome() 
```
- Function to return a valid GNOME string required to create the population. 

```cpp
int cal_fitness(string gnome) 
```
- Function to return the fitness value of a gnome. 
- The fitness value is the path length of the path represented by the GNOME. 

```cpp
int cooldown(int temp) 
```
- Function to return the updated value of the coolong variable


```cpp
bool operator < ( const individual& item )const
```
- Comparator for GNome struct that compares fitness values of 2 chromosomes.
- Result used to decide if mutation is successful.

### TSPUtil
```cpp
void TrojanMap::TSPUtil(std::vector<std::vector<double>> map, double &il) 
```
- Core part of the genetic algorithm.
- Function is set with population size, cooling variable and no. of iterations as parameters
- A while loop is invoked to run until either the cooling variable indicates or the iterations are exhausted.
- Every loop involves creating a new pool of chromosomes / paths based on parent generation.
- Mutation happens by randomly swapping any two genes in the chromosome.
- Once a mutation is established, its fitness is determined and a decision is made whether or not it is a part of the next generation.
- Rejected chromosomes are accepted based on a probabilistic metric.


### TravellingTrojan_GeneticAlgo
```cpp
std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan_GeneticAlgo(std::vector<std::string> &location_ids) 
```
- Function generates the weight matrix based on inpt location ids.
- Invokes TSPUtil.

- **Time Complexity** O(g(nm + nm + n)) with g the number of generations, n the population size and m the size of the individuals
---


### Results
#### Shortest Path Highlighted
##### Target to Popeyes Louisiana Kitchen - A distance of 1.53852 miles
<p align="center"><img src="img/GMaps_target_popeyes.PNG" alt="TargetToPopeyesGM" width="500"/></p>

<p align="center"><img src="img/SP_target_popeyes.PNG" alt="TargetToPopeyes" width="500"/></p>

##### Ralphs to ChickFil-A - A distance of 0.74044 miles
<p align="center"><img src="img/SP_ralphs_chickfila.PNG" alt="RalphsToChickFilA" width="500"/></p>

<p align="center"><img src="img/GMaps_ralphs_chickfila.PNG" alt="R2CGM" width="500"/></p>

##### Travelling Salesman Animation

<p align="center"><img src="img/TSP.gif" alt="TargetToPopeyes" width="500"/></p>

###### Timing Comparison of the different algorithms
- The time taken on our computers by each of the 3 algorithms are given below
 - Travelling salesman using Brute Force - 13.14 seconds
 - Travelling salesman using 2-Opt - 1.435 seconds
 - Travelling salesman using Genetic Algorithm - 4.385 seconds
- By taking the ratios between the 3 algorithms using Brute Force as a baseline we find
 - 2-Opt is 9.156 times faster than Brute Force
 - Genetic Algorithm is 2.996 times faster than Brute Force
 We find that the Genetic Algorithm is faster than the Brute Force. However, it is more accurate than 2-Opt.


---
### Testing the Algorithms
- The algorithms were tested using the tests provided by the TAs for validation of the solution.
- For the heuristic algorithms i.e., 2 opt and Genetic Algorithm, the test cases were established in such way that it passes if an improvement is made or the result remains the same in certain cases.
- Designed tests ourselves to verify our code’s error handling. This was done in 2 ways.
 - Modifying the test cases generated by the TAs.
 - Writing our own test cases to better understand our code.
- Designed test cases to handle corner cases so that outlying inputs/conditions would not affect our final result.
 - For example, tested for Node ID not available in the .CSV file to alert the user.
 - Empty result if input is an invalid location while calculating Shortest Path.
 - Return string as result if input is of size 1.
 - Ensuring that output was the same even when interchanging the destination and source for Shortest Path.

---

### Learnings - _Ashwin Sivakumar_
- Intuition behind shortest path algorithms. Knowledge gained will help take an efficient and informed approach before tackling a problem.
- Power of 2-opt.A good and quicker approximate might come in handy. Not making perfect the enemy of the good!
- Obtained an understanding of the Genetic Algorithm. In addition to understanding how it could be used to solve TSP, I am fascinated by its application in data driven control engineering.
- Handling, organizing and accessing data efficiently and effectively by leveraging the power of STLs.
- Fine tuned my thought process while addressing a problem. Has become more structured and systematic.

### Learnings - _Mackton Vishal Mendonca_

<p align="center"><img src="img/Learnings.gif" alt="Amoeba" width="200"/></p>

- Research on this problem revealed some non computer science related problems that are solved using TSP. 
- (**Physarum polycephalum**, an amoeboid that adapts its morphology to create an efficient path between its food sources).
- Substituting cities and travelling distances for chip components and connecting circuits, the potential of TSP in chip design to integrate more and more transistors on a chip.
- Deeper understanding of object oriented programming for data extraction and usage of structures in C++.
- Understanding the importance of testing an application to not just provide a solution but also catch “not solutions”. To not only catch bugs but also to improve the quality of code.
- This project also helped me understand Unix commands as well as the usage of Github.