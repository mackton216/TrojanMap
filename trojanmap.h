#ifndef TROJAN_MAP_H
#define TROJAN_MAP_H
#define DOT_SIZE 5
#define LINE_WIDTH 3

#include <iostream>
#include <map>
#include <vector>
#include <unordered_set>
#include <set>
#include <algorithm>



// A Node is the location of one point in the map.
class Node {
  public:
    Node(){};
    Node(const Node &n){id = n.id; lat = n.lat; lon = n.lon; name = n.name; neighbors = n.neighbors;};
    std::string id;    // A unique id assign to each point
    double lat;        // Latitude
    double lon;        // Longitude
    std::string name;  // Name of the location. E.g. "Bank of America".
    std::vector<std::string>neighbors;  // List of the ids of all neighbor points.
};

struct individual
{
  std::string gnome;
  double fitness;
  bool operator < ( const individual&  item )const
  {
      return fitness < item.fitness;
  }
};




class TrojanMap {
 public:
  //-----------------------------------------------------
  // TODO: You do not and should not change the following functions:

  // Create the menu.
  void PrintMenu();

  // Read in the data
  void CreateGraphFromCSVFile();

  // Visualization
  // Given a location id, plot the point on the map.
  void PlotPoint(std::string id);

  // Given a lat and lon, plot the point on the map.
  void PlotPoint(double lat, double lon);

  // Given a vector of location ids draws the path (connects the points)
  void PlotPath(std::vector<std::string> &location_ids);

  // Given a vector of location ids draws the points on the map (no path).
  void PlotPoints(std::vector<std::string> &location_ids);

  // Create the videos of the progress to get the path
  void CreateAnimation(std::vector<std::vector<std::string>>);

  // Transform the location to the position on the map
  std::pair<double, double> GetPlotLocation(double lat, double lon);
  //-----------------------------------------------------
  // TODO: Implement these functions and create unit tests for them:

  // Get the Latitude of a Node given its id.
  double GetLat(std::string id);

  // Get the Longitude of a Node given its id.
  double GetLon(std::string id);

  // Get the name of a Node given its id.
  std::string GetName(std::string id);

  // Get the neighbor ids of a Node.
  std::vector<std::string> GetNeighborIDs(std::string id);

  //Extract Node given ID

  Node GetNode(std::string id);

  //Extract Node from name
  Node GetNodeFromName(std::string locname);

  // Get the distance between 2 nodes.
  double CalculateDistance(const Node &a, const Node &b);


  // Calculates the total path length for the locations inside the vector.
  double CalculatePathLength(const std::vector<std::string> &path);

  // Returns a vector of names given a partial name.
  std::vector<std::string> Autocomplete(std::string name);

  // Returns lat and long of the given the name.
  std::pair<double, double> GetPosition(std::string name);

  // Given the name of two locations, it should return the **ids** of the nodes
  // on the shortest path.

  std::vector<std::string> CalculateShortestPath(std::string location1_name,
                                                 std::string location2_name);


  std::vector<double> DijkstraPriorityQueue(int source, std::vector<std::vector<double>> weight_, std::map<int, int> &prev);
  // Given a vector of location ids, it should reorder them such that the path
  // that covers all these points has the minimum length.
  // The return value is a pair where the first member is the total_path,
  // and the second member is the reordered vector of points.
  // (Notice that we don't find the optimal answer. You can return an estimated
  // path.)

double TSPHELPER(int start, int cur_node, double cur_cost, std::vector<std::string> path, std::vector<std::vector<double>> weight_, std::map<int, std::string>& rev_index, std::pair<double,std::vector<std::vector<std::string>>>& progress);

std::pair<double, std::vector<std::vector<std::string>>> TravellingTrojan(std::vector<std::string> &location_ids);


std::vector<int> NN(std::vector<std::vector<double>> weight_, double &dist);
double distance(std::vector<int> pt, int i, std::vector<std::vector<double>> weight_);
std::pair<double, std::vector<std::vector<int>>> TWO_OPT(std::vector<std::vector<double>> weight_, std::pair<double, std::vector<std::vector<int>>> result);

std::pair<double, std::vector<std::vector<std::string>>> TravellingTrojan_2opt(
      std::vector<std::string> &location_ids);
  //-----------------------------------------------------

// Genetic Algorithm

// Comparator for GNOME struct. 
// void TSPUtil(std::vector<std::vector<double>> weight_);
std::pair<double, std::vector<std::vector<std::string>>> TravellingTrojan_GeneticAlgo(std::vector<std::string> &location_ids) ;

void TSPUtil(std::vector<std::vector<double>> map, std::pair<double, std::vector<std::vector<std::string>>>&progress, double &il);

 private:

 // Function to return a random number 
// from start and end 
int rand_num(int start, int end) 
{ 
    int r = end - start; 
    int rnum = start + rand() % r; 
    return rnum; 
} 

// Function to check if the character 
// has already occurred in the string 
bool repeat(std::string s, char ch) 
{ 
    for (int i = 0; i < s.size(); i++) { 
        if (s[i] == ch) 
            return true; 
    } 
    return false; 
} 

// Function to return a mutated GNOME 
// Mutated GNOME is a string 
// with a random interchange 
// of two genes to create variation in species 
std::string mutatedGene(std::string gnome, std::vector<std::vector<double>> map ) 
{   
    int V = map.size();
    while (true) { 
        int r = rand_num(1, V); 
        int r1 = rand_num(1, V); 
        if (r1 != r) { 
            char temp = gnome[r]; 
            gnome[r] = gnome[r1]; 
            gnome[r1] = temp; 
            break; 
        } 
    } 
    return gnome; 
} 

std::string create_gnome(std::vector<std::vector<double>> map) 
{ 
    int V = map.size();
    std::string gnome = "0"; 
    while (true) { 
        if (gnome.size() == V) { 
            gnome += gnome[0]; 
            break; 
        } 
        int temp = rand_num(1, V); 
        if (!repeat(gnome, (char)(temp + 48))) 
            gnome += (char)(temp + 48); 
    } 
    return gnome; 
} 

// Function to return the fitness value of a gnome. 
// The fitness value is the path length 
// of the path represented by the GNOME. 
double cal_fitness(std::string gnome, std::vector<std::vector<double>> map) 
{   
    double f = 0; 
    for (int i = 0; i < gnome.size() - 1; i++) { 
        if (map[gnome[i] - 48][gnome[i + 1] - 48] == INT16_MAX) 
            return INT16_MAX; 
        f += map[gnome[i] - 48][gnome[i + 1] - 48]; 
    } 

    return f; 
} 

// Function to return the updated value 
// of the cooling element. 
int cooldown(int temp) 
{ 
    return (90 * temp) / 100; 
} 
  
  // A map of ids to Nodes.
  std::map<std::string, Node> data;
};

#endif