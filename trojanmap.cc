#include "trojanmap.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <algorithm>
#include <fstream>
#include <locale>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <utility>

#include "opencv2/core.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/videoio.hpp"



//-----------------------------------------------------
// TODO (Students): You do not and should not change the following functions:
//-----------------------------------------------------

/**
 * PrintMenu: Create the menu
 * 
 */
void TrojanMap::PrintMenu() {

  std::string menu =
      "**************************************************************\n"
      "* Select the function you want to execute.                    \n"
      "* 1. Autocomplete                                             \n"
      "* 2. Find the position                                        \n"
      "* 3. CalculateShortestPath - Dijkstra                         \n"
      "* 4. Travelling salesman problem - Brute Force                \n"
      "* 5. Travelling salesman problem - TWO OPT                    \n"
      "* 6. Travelling salesman problem - Genetic Algorithm          \n"
      "* 7. Exit                                                     \n"
      "**************************************************************\n";
  std::cout << menu << std::endl;
  std::string input;
  getline(std::cin, input);
  char number = input[0];
  switch (number)
  {
  case '1':
  {
    menu =
        "**************************************************************\n"
        "* 1. Autocomplete                                             \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "Please input a partial location:";
    std::cout << menu;
    getline(std::cin, input);
    auto results = Autocomplete(input);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    if (results.size() != 0) {
      for (auto x : results) std::cout << x << std::endl;
    } else {
      std::cout << "No matched locations." << std::endl;
    }
    menu = "**************************************************************\n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }
  case '2':
  {
    menu =
        "**************************************************************\n"
        "* 2. Find the position                                        \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "Please input a location:";
    std::cout << menu;
    getline(std::cin, input);
    auto results = GetPosition(input);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    if (results.first != -1) {
      std::cout << "Latitude: " << results.first
                << " Longitude: " << results.second << std::endl;
      PlotPoint(results.first, results.second);
    } else {
      std::cout << "No matched locations." << std::endl;
    }
    menu = "**************************************************************\n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }
  case '3':
  {
    menu =
        "**************************************************************\n"
        "* 3. CalculateShortestPath                                            "
        "      \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "Please input the start location:";
    std::cout << menu;
    std::string input1;
    getline(std::cin, input1);
    menu = "Please input the destination:";
    std::cout << menu;
    std::string input2;
    getline(std::cin, input2);
    auto results = CalculateShortestPath(input1, input2);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    if (results.size() != 0) {
      for (auto x : results) std::cout << x << std::endl;
      PlotPath(results);
    } else {
      std::cout << "No route from the start point to the destination."
                << std::endl;
    }
    menu = "**************************************************************\n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }

  case '4':
  {
    menu =
        "**************************************************************\n"
        "* 4. Travelling salesman problem - Brute Force                \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "In this task, we will select N random points on the map and you need to find the path to travel these points and back to the start point.";
    std::cout << menu << std::endl << std::endl;
    menu = "Please input the number of the places:";
    std::cout << menu;
    getline(std::cin, input);
    int num = std::stoi(input);
    std::vector<std::string> keys;
    for (auto x : data) {
      keys.push_back(x.first);
    }
    std::vector<std::string> locations;
    srand(time(NULL));
    for (int i = 0; i < num; i++)
      locations.push_back(keys[rand() % keys.size()]);
    PlotPoints(locations);
    std::cout << "Calculating ..." << std::endl;
    auto results = TravellingTrojan(locations);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    CreateAnimation(results.second);
    if (results.second.size() != 0) {
      for (auto x : results.second[results.second.size()-1]) std::cout << x << std::endl;
      menu = "**************************************************************\n";
      std::cout << menu;
      std::cout << "The distance of the path is:" << results.first << std::endl;
      PlotPath(results.second[results.second.size()-1]);
    } else {
      std::cout << "The size of the path is 0" << std::endl;
    }
    menu = "**************************************************************\n"
           "You could find your animation at src/lib/output.avi.          \n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }

  case '5':
  {
    menu =
        "**************************************************************\n"
        "* 5. Travelling salesman problem  - TWO OPT                   \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "In this task, we will select N random points on the map and you need to find the path to travel these points and back to the start point using 2-opt.";
    std::cout << menu << std::endl << std::endl;
    menu = "Please input the number of the places:";
    std::cout << menu;
    getline(std::cin, input);
    int num = std::stoi(input);
    std::vector<std::string> keys;
    for (auto x : data) {
      keys.push_back(x.first);
    }
    std::vector<std::string> locations;
    srand(time(NULL));
    for (int i = 0; i < num; i++)
      locations.push_back(keys[rand() % keys.size()]);
    PlotPoints(locations);
    std::cout << "Calculating ..." << std::endl;
    auto results = TravellingTrojan_2opt(locations);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    CreateAnimation(results.second);
    if (results.second.size() != 0) {
      for (auto x : results.second[results.second.size()-1]) std::cout << x << std::endl;
      menu = "**************************************************************\n";
      std::cout << menu;
      std::cout << "The distance of the path is:" << results.first << std::endl;
      PlotPath(results.second[results.second.size()-1]);
    } else {
      std::cout << "The size of the path is 0" << std::endl;
    }
    menu = "**************************************************************\n"
           "You could find your animation at src/lib/output.avi.          \n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }

  case '6':
  {
    menu =
        "**************************************************************\n"
        "* 6. Travelling salesman problem  - Genetic Algorithm         \n"
        "**************************************************************\n";
    std::cout << menu << std::endl;
    menu = "In this task, we will select N random points on the map and you need to find the path to travel these points and back to the start point using 2-opt.";
    std::cout << menu << std::endl << std::endl;
    menu = "Please input the number of the places:";
    std::cout << menu;
    getline(std::cin, input);
    int num = std::stoi(input);
    std::vector<std::string> keys;
    for (auto x : data) {
      keys.push_back(x.first);
    }
    std::vector<std::string> locations;
    srand(time(NULL));
    for (int i = 0; i < num; i++)
      locations.push_back(keys[rand() % keys.size()]);
    PlotPoints(locations);
    std::cout << "Calculating ..." << std::endl;
    auto results = TravellingTrojan_GeneticAlgo(locations);
    menu = "*************************Results******************************\n";
    std::cout << menu;
    CreateAnimation(results.second);
    if (results.second.size() != 0) 
    {
      // for (auto x : results.second[results.second.size()-1]) std::cout << x << std::endl;
      menu = "**************************************************************\n";
      std::cout << menu;
      std::cout << "The distance of the path is:" << results.first << std::endl;
      PlotPath(results.second[results.second.size()-1]);
    }
     else
      {
      std::cout << "The size of the path is 0" << std::endl;
      }
    menu = "**************************************************************\n"
           "You could find your animation at src/lib/output.avi.          \n";
    std::cout << menu << std::endl;
    PrintMenu();
    break;
  }

  case '7':
    break;
  default:
    std::cout << "Please select 1 - 5." << std::endl;
    PrintMenu();
    break;
  }
}


/**
 * CreateGraphFromCSVFile: Read the map data from the csv file
 * 
 */
void TrojanMap::CreateGraphFromCSVFile() {
  std::fstream fin;
  fin.open("src/lib/map.csv", std::ios::in);
  std::string line, word;

  getline(fin, line);
  while (getline(fin, line)) {
    std::stringstream s(line);

    Node n;
    int count = 0;
    while (getline(s, word, ',')) {
      word.erase(std::remove(word.begin(), word.end(), '\''), word.end());
      word.erase(std::remove(word.begin(), word.end(), '"'), word.end());
      word.erase(std::remove(word.begin(), word.end(), '['), word.end());
      word.erase(std::remove(word.begin(), word.end(), ']'), word.end());
      if (count == 0)
        n.id = word;
      else if (count == 1)
        n.lat = stod(word);
      else if (count == 2)
        n.lon = stod(word);
      else if (count == 3)
        n.name = word;
      else {
        word.erase(std::remove(word.begin(), word.end(), ' '), word.end());
        n.neighbors.push_back(word);
      }
      count++;
    }
    data[n.id] = n;
  }
  fin.close();
}

/**
 * PlotPoint: Given a location id, plot the point on the map
 * 
 * @param  {std::string} id : location id
 */
void TrojanMap::PlotPoint(std::string id) {
  std::string image_path = cv::samples::findFile("src/lib/input.jpg");
  cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
  auto result = GetPlotLocation(data[id].lat, data[id].lon);
  cv::circle(img, cv::Point(result.first, result.second), DOT_SIZE,
             cv::Scalar(0, 0, 255), cv::FILLED);
  cv::imshow("TrojanMap", img);
  cv::waitKey(1);
}
/**
 * PlotPoint: Given a lat and a lon, plot the point on the map
 * 
 * @param  {double} lat : latitude
 * @param  {double} lon : longitude
 */
void TrojanMap::PlotPoint(double lat, double lon) {
  std::string image_path = cv::samples::findFile("src/lib/input.jpg");
  cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
  auto result = GetPlotLocation(lat, lon);
  cv::circle(img, cv::Point(int(result.first), int(result.second)), DOT_SIZE,
             cv::Scalar(0, 0, 255), cv::FILLED);
  cv::startWindowThread();
  cv::imshow("TrojanMap", img);
  cv::waitKey(1);
}

/**
 * PlotPath: Given a vector of location ids draws the path (connects the points)
 * 
 * @param  {std::vector<std::string>} location_ids : path
 */
void TrojanMap::PlotPath(std::vector<std::string> &location_ids) {
  std::string image_path = cv::samples::findFile("src/lib/input.jpg");
  cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
  auto start = GetPlotLocation(data[location_ids[0]].lat, data[location_ids[0]].lon);
  cv::circle(img, cv::Point(int(start.first), int(start.second)), DOT_SIZE,
             cv::Scalar(0, 0, 255), cv::FILLED);
  for (auto i = 1; i < location_ids.size(); i++) {
    auto start = GetPlotLocation(data[location_ids[i - 1]].lat, data[location_ids[i - 1]].lon);
    auto end = GetPlotLocation(data[location_ids[i]].lat, data[location_ids[i]].lon);
    cv::circle(img, cv::Point(int(end.first), int(end.second)), DOT_SIZE,
               cv::Scalar(0, 0, 255), cv::FILLED);
    cv::line(img, cv::Point(int(start.first), int(start.second)),
             cv::Point(int(end.first), int(end.second)), cv::Scalar(0, 255, 0),
             LINE_WIDTH);
  }
  cv::startWindowThread();
  cv::imshow("TrojanMap", img);
  cv::waitKey(1);
}

/**
 * PlotPoints: Given a vector of location ids draws the points on the map (no path).
 * 
 * @param  {std::vector<std::string>} location_ids : points
 */
void TrojanMap::PlotPoints(std::vector<std::string> &location_ids) {
  std::string image_path = cv::samples::findFile("src/lib/input.jpg");
  cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
  for (auto x : location_ids) {
    auto result = GetPlotLocation(data[x].lat, data[x].lon);
    cv::circle(img, cv::Point(result.first, result.second), DOT_SIZE,
               cv::Scalar(0, 0, 255), cv::FILLED);
  }
  cv::imshow("TrojanMap", img);
  cv::waitKey(1);
}


/**
 * CreateAnimation: Create the videos of the progress to get the path
 * 
 * @param  {std::vector<std::vector<std::string>>} path_progress : the progress to get the path
 */
void TrojanMap::CreateAnimation(std::vector<std::vector<std::string>> path_progress){
  cv::VideoWriter video("src/lib/output.avi", cv::VideoWriter::fourcc('M','J','P','G'), 10, cv::Size(1248,992));
  for(auto location_ids: path_progress) {
    std::string image_path = cv::samples::findFile("src/lib/input.jpg");
    cv::Mat img = cv::imread(image_path, cv::IMREAD_COLOR);
    auto start = GetPlotLocation(data[location_ids[0]].lat, data[location_ids[0]].lon);
    cv::circle(img, cv::Point(int(start.first), int(start.second)), DOT_SIZE,
              cv::Scalar(0, 0, 255), cv::FILLED);
    for (auto i = 1; i < location_ids.size(); i++) {
      auto start = GetPlotLocation(data[location_ids[i - 1]].lat, data[location_ids[i - 1]].lon);
      auto end = GetPlotLocation(data[location_ids[i]].lat, data[location_ids[i]].lon);
      cv::circle(img, cv::Point(int(end.first), int(end.second)), DOT_SIZE,
                cv::Scalar(0, 0, 255), cv::FILLED);
      cv::line(img, cv::Point(int(start.first), int(start.second)),
              cv::Point(int(end.first), int(end.second)), cv::Scalar(0, 255, 0),
              LINE_WIDTH);
    }
    video.write(img);
    cv::imshow("TrojanMap", img);
    cv::waitKey(1);
  }
	video.release();
}
/**
 * GetPlotLocation: Transform the location to the position on the map
 * 
 * @param  {double} lat         : latitude 
 * @param  {double} lon         : longitude
 * @return {std::pair<double, double>}  : position on the map
 */
std::pair<double, double> TrojanMap::GetPlotLocation(double lat, double lon) {
  std::pair<double, double> bottomLeft(34.01001, -118.30000);
  std::pair<double, double> upperRight(34.03302, -118.26502);
  double h = upperRight.first - bottomLeft.first;
  double w = upperRight.second - bottomLeft.second;
  std::pair<double, double> result((lon - bottomLeft.second) / w * 1248,
                                   (1 - (lat - bottomLeft.first) / h) * 992);
  return result;
}

//-----------------------------------------------------
// TODO: Student should implement the following:
//-----------------------------------------------------
/**
 * GetLat: Get the latitude of a Node given its id.
 * 
 * @param  {std::string} id : location id
 * @return {double}         : latitude
 */
double TrojanMap::GetLat(std::string id) 
{ 
 return data[id].lat;
}

/**
 * GetLon: Get the longitude of a Node given its id. 
 * 
 * @param  {std::string} id : location id
 * @return {double}         : longitude
 */


double TrojanMap::GetLon(std::string id) 
{  
   return data[id].lon;
}

/**
 * GetName: Get the name of a Node given its id.
 * 
 * @param  {std::string} id : location id
 * @return {std::string}    : name
 */
std::string TrojanMap::GetName(std::string id) 
{ 
  return data[id].name;
}

/**
 * GetNeighborIDs: Get the neighbor ids of a Node.
 * 
 * @param  {std::string} id            : location id
 * @return {std::vector<std::string>}  : neighbor ids
 */
std::vector<std::string> TrojanMap::GetNeighborIDs(std::string id) 
{   

    std::vector<std::string> result;

    for(int i = 0; i < data[id].neighbors.size(); i++)
    {
      result.push_back(data[id].neighbors[i]);
    }
  return result;
} 


/**
 * CalculateDistance: Get the distance between 2 nodes. 
 * 
 * @param  {Node} a  : node a
 * @param  {Node} b  : node b
 * @return {double}  : distance in mile
 */
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
  
  // TODO: Use Haversine Formula:
  // dlon = lon2 - lon1;
  // dlat = lat2 - lat1;
  // a = (sin(dlat / 2)) ^ 2 + cos(lat1) * cos(lat2) * (sin(dlon / 2)) ^ 2;
  // c = 2 * arcsin(min(1, sqrt(a)));
  // distances = 3961 * c;

  // where 3961 is the approximate radius of the earth at the latitude of
  // Washington, D.C., in miles
}

/**
 * CalculatePathLength: Calculates the total path length for the locations inside the vector.
 * 
 * @param  {std::vector<std::string>} path : path
 * @return {double}                        : path length
 */

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

/**
 * Autocomplete: Given a parital name return all the possible locations with
 * partial name as the prefix.
 *
 * @param  {std::string} name          : partial name
 * @return {std::vector<std::string>}  : a vector of full names
 */

std::vector<std::string> TrojanMap::Autocomplete(std::string name) {

  int n = name.size();
  std::string temp;
  std::vector<std::string> results = {};
  
  std::transform(name.begin(), name.end(), name.begin(),
    [](unsigned char c){ return std::tolower(c); });

  std::map<std::string, Node>::iterator it;

  for(it = data.begin(); it!=data.end();it++)
  { 

      if(it->second.name.empty() == 0)
      {
           
        temp.append(it->second.name.begin(),it->second.name.begin() + n);

        std::transform(temp.begin(), temp.end(), temp.begin(),
        [](unsigned char c1){ return std::tolower(c1); });

        if(name.compare(temp)==0)
          {
            results.push_back(it->second.name);
          }
      }
      temp.clear();
  }
  return results;
}

/**
 * GetPosition: Given a location name, return the position.
 *
 * @param  {std::string} name          : location name
 * @return {std::pair<double,double>}  : (lat, lon)
 */
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


std::vector<double> TrojanMap::DijkstraPriorityQueue(int source, std::vector<std::vector<double>> weight_, std::map<int, int> &prev)
{
  std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>>q;
  std::vector<double> d(weight_.size(), DBL_MAX), d_temp(weight_.size(), DBL_MAX);

  d[source] = 0;
  q.push(std::make_pair(0, source));

  while(!q.empty())
  { 

    int u = q.top().second;
    q.pop();
    for(int j = 0; j < weight_.size(); j++)
    {
      if(d[j] > d[u] + weight_[u][j])
      {
        d[j] = d[u] + weight_[u][j];
        q.push(std::make_pair(d[j], j));
        prev[j] = u;
      }
    }
  }
  return d;
}


/**
 * CalculateShortestPath: Given 2 locations, return the shortest path which is a
 * list of id.
 *
 * @param  {std::string} location1_name     : start
 * @param  {std::string} location2_name     : goal
 * @return {std::vector<std::string>}       : path
 * 
 */


std::vector<std::string> TrojanMap::CalculateShortestPath(
  
  std::string location1_name, std::string location2_name) {
  std::vector<std::string> path = {};
  Node a, b; 
  double dist;
  int m = 0;
  std::map<std::string, Node>::iterator it;

  int flag = 0;
  for(it = data.begin(); it!=data.end();it++)
  {
    if(it->second.name.compare(location2_name) == 0)
    {
      flag = 1;
      break;         
    }
  }

  if(flag == 0)
  {
    path = {};
    return path;
  }


  a = GetNodeFromName(location1_name);
  b = GetNodeFromName(location2_name); 
  
  std::vector<std::vector<double>> weight_(2237, std::vector<double>(2237, DBL_MAX));

  std::map<std::string, Node>::iterator itr1, itr2;

  std::map<std::string, int>index;
  std::map<std::string, int>::iterator ie;
  std::map<int, std::string> rev_in;

  

  std::map<int, int> prev;

  for(itr1 = data.begin(); itr1!=data.end(); itr1++) // current node
  {
    index[itr1->first] = m;
    m++;
  }

   for(ie = index.begin(); ie != index.end(); ie++)
  {
      rev_in[ie->second] = ie->first;
  }

  int source = index[a.id];
  int destination = index[b.id];

  if(source == destination)
  {
    path.push_back(rev_in[source]);
    return path;
  }

 
  int found = 0, mx = 0;

  for(itr1 = data.begin(); itr1!=data.end(); itr1++) // current node
  { 

    for(itr2= data.begin(); itr2!=data.end(); itr2++) // iterate through all the nodes
    {
      for(int k = 0; k < itr1->second.neighbors.size(); k++) // check neighbors of current node
      {
        if(found == 1)
        {
          continue;
        }
        if(itr1->second.neighbors[k].compare(itr2->first) == 0) // check if neighborID and node ID match, if so find distance and update AM
        {
          found = 1;
          dist = CalculateDistance(itr1->second, itr2->second);
          weight_[index[itr1->first]][index[itr2->first]] = dist;
        }
      }
      if(found == 0)
      {
        weight_[index[itr1->first]][index[itr2->first]] = DBL_MAX;
      }
      found = 0;
    }

  }


   for(itr1 = data.begin(); itr1!=data.end(); itr1++) // diagonal elements set to zero
  {
    weight_[index[itr1->first]][index[itr1->first]] = 0;
  }

  for(int i = 0; i < 2237; i++)
  {
      if(weight_[source][i]!= DBL_MAX)
      {
        prev[i] = source;
      }
      else
      {
        prev[i] = DBL_MAX;
      }
  }
 
  std::vector<double> x = DijkstraPriorityQueue(source, weight_, prev);  

    
  int mi = destination;

  path.push_back(rev_in[destination]);
  while(mi!=source)
  { 
    path.push_back(rev_in[prev[mi]]);
    mi = prev[mi];
  }

  std::reverse(path.begin(),path.end());
  return path;
  
}


/**
 * Travelling salesman problem: Given a list of locations, return the shortest
 * path which visit all the places and back to the start point.
 *
 * @param  {std::vector<std::string>} input : a list of locations needs to visit
 * @return {std::pair<double, std::vector<std::vector<std::string>>} : a pair of total distance and the all the progress to get final path
 */


double TrojanMap::TSPHELPER(int start, int cur_node, double cur_cost, std::vector<std::string> path, std::vector<std::vector<double>> weight_, std::map<int, std::string>& rev_index, std::pair<double, std::vector<std::vector<std::string>>>&progress)
{
 
  double result = INT16_MAX;

    path.push_back(rev_index[cur_node]);

    if(path.size() == weight_.size())
    {
      path.push_back(rev_index[start]);
      progress.second.push_back(path);
      return cur_cost + weight_[cur_node][start];
    }

    for(int i = 0; i < weight_.size(); i++)
    {
      if(i!=cur_node && std::find(path.begin(), path.end(), rev_index[i]) == path.end())
      { 
        result = std::min(result, TSPHELPER(start, i, cur_cost + weight_[cur_node][i], path, weight_, rev_index, progress));
      }
    }

  progress.first = result;
  return result;
}




std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan(std::vector<std::string> &location_ids) 
{
 
  std::pair<double, std::vector<std::vector<std::string>>> progress;
  std::map<std::string, int>index;
  std::map<int, std::string>rev_index;

  
  double dist;
  std::vector<std::vector<double>> weight_(location_ids.size(), std::vector<double> (location_ids.size()));
  // create indices

  for(int i = 0; i < location_ids.size(); i++)// current node
  {
    index[location_ids[i]] = i;
  }

  for(int i = 0; i < location_ids.size(); i++)// current node
  {
    rev_index[i] = location_ids[i];
  }
  
  // generate weight matrices
  for(int i = 0; i < location_ids.size(); i++)
  { 
    for(int j = 0; j < location_ids.size(); j++)
    {
      if(i == j)
      {
        weight_[i][j] = 0;
      }

      else
      {      
      dist = CalculateDistance(data[location_ids[i]], data[location_ids[j]]);      
      weight_[i][j] = dist;
      }
    }
  }
  
  std::vector<std::string> path;

  progress.first = 0;
  int start = 0;
  double d;
  double rs = TSPHELPER(start, start, 0, path, weight_, rev_index, progress);
  progress.first = rs;

  for(int i = 0; i < progress.second.size(); i++)
  {
    d = CalculatePathLength(progress.second[i]);
    if(d == progress.first)
    {
    progress.second.push_back(progress.second[i]);
    progress.second.erase(progress.second.begin()+i);
    }
  }

  
  return progress;
} 


std::vector<int> TrojanMap::NN(std::vector<std::vector<double>> weight_, double &dist)
{
  std::vector<int> visited;
  std::vector<int> path;
  int cur = 0, j, start = 0;

  path.push_back(cur);
  visited.push_back(cur);
  
  
  int min_id;
  double min;  

 
  while(visited.size() < weight_.size())
  {
    min = DBL_MAX;
    for(j = 0; j < weight_.size(); j++)
    {
      if(std::find(visited.begin(), visited.end(), j) == visited.end())
      {
        if(weight_[cur][j] < min)
        {
          min = weight_[cur][j];
          min_id = j;
        }
      }
    }
    dist += weight_[cur][min_id];
    cur = min_id;
    path.push_back(min_id);
    visited.push_back(min_id);
  }
  dist+=weight_[cur][start];

  path.push_back(start);

  return path;
}


double TrojanMap::distance(std::vector<int> pt, int i, std::vector<std::vector<double>> weight_)
{
  double dist = 0;
  int start = i;

  for(int j = 0; j < weight_.size() - 1;j++)
  {
    dist += weight_[pt[i]][pt[i+1]];
    i++;
  }

  dist+=weight_[pt[i]][pt[start]];

 
  return dist;
}


std::pair<double, std::vector<std::vector<int>>> TrojanMap::TWO_OPT(std::vector<std::vector<double>> weight_, std::pair<double, std::vector<std::vector<int>>> result)
{
  
  std::vector<int> shortest = result.second[0];
  std::vector<int> s_temp;
  std::vector<int> mp, m_temp;
  std::vector<int> temp_store;
  int j = 0, p;
  double dist = result.first;

  int flag = 1,x=0;
  while(flag == 1)
  {
    mp = result.second[result.second.size()-1];
    mp.pop_back();
    mp.insert(std::end(mp), std::begin(mp), std::end(mp));
    s_temp = shortest;
    flag = 0;
    for(j = 0; j < weight_.size(); j++)
    { 
      m_temp = mp;
      std::reverse(m_temp.begin() + j + 1, m_temp.begin() + j + 4);

      // calculate new distance
      dist = distance(m_temp, j, weight_);
      if(dist < result.first)
      {
        flag = 1;
        
        result.first = dist;
        p = j;
        for(int k = 0; k < weight_.size(); k++)
        {
          temp_store.push_back(m_temp[p]);
          p++;
        }
        temp_store.push_back(m_temp[j]);
        result.second.push_back(temp_store);
        shortest = temp_store;        
        temp_store.clear();
        break;
      }
    }
  }

return result;
}

std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan_2opt(std::vector<std::string> &location_ids)
{
  std::pair<double, std::vector<std::vector<std::string>>> progress;
  std::map<std::string, int>index;
  std::map<int, std::string>rev_index;
  

  std::vector<std::string> loc_temp = location_ids;
  loc_temp.push_back(location_ids[0]);

  // double d_check = CalculatePathLength(loc_temp);
  // std::cout<<d_check<<"\n";
  

  double d = 0;  

  std::vector<std::vector<double>> weight_(location_ids.size(), std::vector<double> (location_ids.size()));

  for(int i = 0; i < location_ids.size(); i++)// current node
  {
    index[location_ids[i]] = i;
  }


  for(int i = 0; i < location_ids.size(); i++)// current node
  {
    rev_index[i] = location_ids[i];
  }
  
  // generate weight matrices
  for(int i = 0; i < location_ids.size(); i++)
  { 
    for(int j = 0; j < location_ids.size(); j++)
    {
      if(i == j)
      {
        weight_[i][j] = 0;
      }

      else
      {      
      d = CalculateDistance(data[location_ids[i]], data[location_ids[j]]);      
      weight_[i][j] = d;
      }
    }
  }

  // Start 2-opt


  double path_dist = 0;
  
  //Obtain nearest neighbor heuristic before starting 2 opt

  auto nn = NN(weight_, path_dist);
  std::cout<<"\nDistance before 2 opt"<<" "<<path_dist<<"miles\n";
  std::pair<double, std::vector<std::vector<int>>>result;
 

  std::vector<int> path_temp;

  result.first = path_dist;
  result.second.push_back(nn);
  

  // Call 2opt 

  result = TWO_OPT(weight_, result);


  int z = result.second.size();
  int o = weight_.size() + 1;
  std::vector<std::vector<std::string>> path_history;

  std::vector<std::string> temp_path;
  std::vector<int> temp_ind;

  for(int u = 0; u < result.second.size(); u++)
  {
    temp_ind = result.second[u];
    
    for(int l = 0; l < temp_ind.size(); l++)
    {
      temp_path.push_back(rev_index[temp_ind[l]]);
    }
    path_history.push_back(temp_path);
    temp_path.clear();
    temp_ind.clear();
  }


  progress.first = result.first;

  std::cout<<"Distance after 2 opt"<<" "<<progress.first<<" miles\n";
  std::cout<<"Improvement of"<<" "<<path_dist - progress.first<<" miles\n\n";
  progress.second = path_history;
  return progress;
}







// Genetic Algorithm
void TrojanMap::TSPUtil(std::vector<std::vector<double>> map, std::pair<double, std::vector<std::vector<std::string>>> &progress, double &il) 
{ 
    // Generation Number 
    int gen = 1; 
    // Number of Gene Iterations 
    int gen_thres = 15; 
    // Population Size
    int POP_SIZE = 10;
    
    
    std::vector<struct individual> population; 
    struct individual temp; 
  
    // Populating the GNOME pool. 
    for (int i = 0; i < POP_SIZE; i++) { 
        temp.gnome = create_gnome(map); 
        temp.fitness = cal_fitness(temp.gnome, map); 
        population.push_back(temp); 
    } 
  
    std::cout << "\nInitial population: " << std::endl 
         << "GNOME     FITNESS VALUE\n"; 
  
    for (int i = 0; i < POP_SIZE; i++)
    { 
        std::cout << population[i].gnome << "   "
             << population[i].fitness << std::endl;
        il+=population[i].fitness; 
    } 
    il = il / population.size();
    std::cout << "\n"; 
  
    bool found = false; 
    int temperature = 10000; 
  
    // Iteration to perform 
    // population crossing and gene mutation. 
    while (temperature > 1000 && gen <= gen_thres) { 
        std ::sort(population.begin() , population.end() ,[](individual ptr_l , individual ptr_r) { return ptr_l < ptr_r;} );
        std::cout << "\nCurrent temp: " << temperature << "\n"; 
        std::vector<struct individual> new_population; 
  
        for (int i = 0; i < POP_SIZE; i++) { 
            struct individual p1 = population[i]; 
          
            while (true) { 
               
                std::string new_g = mutatedGene(p1.gnome, map); 
                struct individual new_gnome; 
                new_gnome.gnome = new_g; 
                new_gnome.fitness = cal_fitness(new_gnome.gnome, map); 
                
                if (new_gnome.fitness <= population[i].fitness) { 
                    new_population.push_back(new_gnome); 
                    break; 
                } 
                else { 
                      // Accepting the rejected children at 
                    // a possible probablity above threshold. 
                    float prob = pow(2.7, 
                                     -1 * ((float)(new_gnome.fitness 
                                                   - population[i].fitness) 
                                           / temperature)); 
                    if (prob > 0.5) { 
                        new_population.push_back(new_gnome); 
                        break; 
                    } 
                } 
            } 
        } 
  
        temperature = cooldown(temperature); 
        population = new_population; 
        std::cout << "Generation " << gen << " \n"; 
        std::cout << "GNOME     FITNESS VALUE\n"; 
        double min = population[0].fitness;
        int min_gnome;
        for (int i = 0; i < POP_SIZE; i++)
        { 
          std::cout << population[i].gnome << "   "<< population[i].fitness << std::endl; 
          if(population[i].fitness < min)
          {
            min = population[i].fitness;
            min_gnome = i;
          }
        }
        std::vector<std::string> x = {};
        x.push_back(population[min_gnome].gnome);
        if(min < progress.first)
        {
          progress.first = min;
        }
        progress.second.push_back(x);

        gen++; 
    } 
} 



std::pair<double, std::vector<std::vector<std::string>>> TrojanMap::TravellingTrojan_GeneticAlgo(std::vector<std::string> &location_ids) 

{
 
  std::pair<double, std::vector<std::vector<std::string>>> progress;
  std::pair<double, std::vector<std::vector<std::string>>> result;
  progress.first = DBL_MAX;
  std::map<std::string, int>index;
  std::map<int, std::string>rev_index;

  
  double dist;
  std::vector<std::vector<double>> weight_(location_ids.size(), std::vector<double> (location_ids.size()));
  // create indices

  for(int i = 0; i < location_ids.size(); i++)// current node
  {
    index[location_ids[i]] = i;
  }

  for(int i = 0; i < location_ids.size(); i++)// current node
  {
    rev_index[i] = location_ids[i];
  }
  
  // generate weight matrices
  for(int i = 0; i < location_ids.size(); i++)
  { 
    for(int j = 0; j < location_ids.size(); j++)
    {
      if(i == j)
      {
        weight_[i][j] = 0;
      }

      else
      {      
      dist = CalculateDistance(data[location_ids[i]], data[location_ids[j]]);      
      weight_[i][j] = dist;
      }
    }
  }
  double il = 0;
  TSPUtil(weight_, progress, il);

  result.first = progress.first;
  std::cout<<"Initial population average"<<" "<<il<<" miles\n";
  std::cout<<"Distance after multiple mutations"<<" "<<result.first<<" miles\n";
  std::cout<<"Improvement of"<<" "<<il-result.first<<" miles\n\n";

  
  for(int i = 0; i < progress.second.size(); i++)
  {
    std::vector<std::string> x;
    for(int j = 0; j < progress.second[i].size(); j++)
    {
      for(int k = 0; k < progress.second[i][j].size(); k++)
      {
        int z = (int)progress.second[i][j].at(k) - 48;
        x.push_back(rev_index[z]);
      }
    }
    result.second.push_back(x);
  }

  return result;
}