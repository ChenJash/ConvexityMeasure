#include <iostream>
#include <vector>
#include <ctime>
#include "../convexity/measureMoS.h"

using namespace std;

// g++ test_mos.cpp -std=c++11

void draw(std::vector<int> & grid, int label, int x1, int y1, int x2, int y2, int width){
    for(auto i = x1; i <= x2; i++){
        for(auto j = y1; j <= y2; j++){
            int temp = i * width + j;
            grid[temp] = label;
        }
    }
}

int main() {
    //// test polygon 
    // PointList polygon;
    // polygon.push_back({0,0});
    // polygon.push_back({1,0});
    // // polygon.push_back({1,0});
    // polygon.push_back({1,1});
    // polygon.push_back({0,0.1});
    // polygon.push_back({-1,1});
    // polygon.push_back({-1,0});
    // cout << checkPolygonConvexByMoS(polygon, vector<PointList>())[0] << endl;
    // return 0;
    
    //// test grid layout
    std::vector<int> grid_asses;
    std::vector<int> cluster_labels;
    for(auto i = 0; i < 400; i++) {
        grid_asses.push_back(i);
    }
    for(auto i = 0; i < 400; i++) {
        cluster_labels.push_back(0);
    }
    draw(cluster_labels, 1, 1, 1, 3, 9, 20);

    auto start = clock();
    auto res = checkConvexForMoS(grid_asses, cluster_labels);
    cout << res[0] / res[1] << endl << "Time:" << clock() - start << "ms" <<endl;
}