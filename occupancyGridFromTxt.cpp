// i need to learn to read from standard input in c++
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

typedef std::vector< std::vector<int>> occupancyGridType; 
occupancyGridType makeOccupancyGrid(int x_dim, int y_dim, std::vector<std::string> rowVector);

int main() {
    // need variables all definied here at the top
    std::string line;
    std::ifstream myfile("../test_environments/grid_envs/environment1.txt");
    int x_dim, y_dim; 
    float x_start, y_start, x_goal, y_goal;
    
    // vector for each string
    std::vector<std::string> rowVector; 
    std::string currentRow;

    if (myfile.is_open()) {
        int i = 0;
        while (std::getline(myfile, line)) {
            std::stringstream streamof(line);
            if (i == 0){
                std::cout << line << std::endl;
                streamof >> x_dim;  // lets worry about htis later
            }
            else if (i==1){
                streamof >> y_dim;
            }
            else if ((i>1) && (i<=11)) {  // hard coding in a 10x10 
                //streamof >> currentRow; 
                std::cout << line << std::endl;
                rowVector.push_back(line);
                // i probably process this into like a 10x10 array right?
                // i think so
            }
            else if (i==12) {
                streamof>> x_start; 
            }
            else if (i==13) {
                streamof>> y_start;
            }
            else if (i==14) {
                streamof >> x_goal;
            }
            else if (i==15) {
                streamof >> y_goal;
            }
            else {
                std::cout<< "finished processing the input file" << std::endl;
            }
            i++;
        }
        myfile.close();
    } 
    std::cout << "the great results were that: " << std::endl;
    std::cout << "x_dim = " <<x_dim << std::endl;
    std::cout << "y_dim = " <<y_dim << std::endl;
    std::cout << "x_start= " <<x_start<< std::endl;
    std::cout << "x_start= " <<x_start<< std::endl;
    std::cout << "y_start= " <<y_start<< std::endl;
    std::cout << "x_goal= " <<x_goal<< std::endl;
    std::cout << "y_goal= " <<y_goal<< std::endl;
    std::cout << "trying to make occupancy grid" << std::endl;
    int grid[y_dim][x_dim];
    occupancyGridType occupancyGrid = makeOccupancyGrid(x_dim, y_dim, rowVector);
    for (auto &r : occupancyGrid) {
        for (auto &z : r) {
            std::cout << z;
        }
        std::cout << std::endl;
    }
    return 0;
};

occupancyGridType makeOccupancyGrid(int x_dim, int y_dim, std::vector<std::string> rowVector){
    // we have a vector of strings, but we want to turn it into a vector of ints
    std::vector<std::vector<int>> occupancyGrid;
    std::vector<int> occupancyGridRow;
    //int grid[y_dim][x_dim];
    for (auto &r : rowVector) { 
        //std::cout << r << std::endl;
        for (auto &t : r) {                     // iterates over vector elements right?
            if (t != ' ') {                      // 32 indicates a space
                //std::cout << t ;
                if (t == '-') {
                    occupancyGridRow.push_back(0);
                    //std::cout << "free space";
                }
                else if (t == '#'){ 
                    occupancyGridRow.push_back(1);
                    //std::cout << "obstacle"; 
                }
            }
        } // iteration over valid row elements 
        occupancyGrid.push_back(occupancyGridRow);
        occupancyGridRow.clear();
    };
    //std::cout << std::endl;
    //std::cout << grid[9][1] << std::endl;
    return occupancyGrid;
};
