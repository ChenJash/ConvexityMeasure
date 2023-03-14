#include <iostream>
#include <vector>
#include <ctime>
#include "../convexity/measureAlphaTriples.h"

using namespace std;

// g++ test_atriple.cpp -std=c++11

void draw(std::vector<int> & grid, int label, int x1, int y1, int x2, int y2, int width){
    for(auto i = x1; i <= x2; i++){
        for(auto j = y1; j <= y2; j++){
            int temp = i * width + j;
            grid[temp] = label;
        }
    }
}

int main() {
    //// test grid layout
    std::vector<int> grid_asses;
    std::vector<int> cluster_labels;
    for(auto i = 0; i < 400; i++) {
        grid_asses.push_back(i);
    }
    for(auto i = 0; i < 400; i++) {
        cluster_labels.push_back(0);
    }
    draw(cluster_labels, 1, 1, 1, 9, 9, 20);

    auto start = clock();
    auto res = checkConvexForAlphaTriples(grid_asses, cluster_labels, 0.25);
    cout << res[0] / res[1] << endl << "Time:" << clock() - start << "ms" <<endl;
}