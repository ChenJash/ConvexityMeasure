#include "../convexity/measureCV.h"
#include <iostream>
#include <vector>
using namespace std;

void test1() {
    int N; cin >> N;
    int num; cin >> num;
    int square_len; cin >> square_len;
    int maxLabel; cin >> maxLabel;

    int *grid_asses = new int[N];
    int *cluster_labels = new int[num];
    for (int i = 0; i < N; i++) {
        cin >> grid_asses[i];
    }
    for (int i = 0; i < num; i++) {
        cin >> cluster_labels[i];
    }

    auto ret = checkConvexForCVArray(grid_asses, cluster_labels, N, num, square_len, maxLabel);
    cout << ret[0] << ", " << ret[1] << endl;

    delete [] grid_asses;
    delete [] cluster_labels;
}

void test2() {
    int N; cin >> N;
    vector<int> _grid_asses(N);
    vector<int> _cluster_labels(N);
    for (int i = 0; i < _grid_asses.size(); i++) {
        cin >> _grid_asses[i];
    }
    for (int i = 0; i < _cluster_labels.size(); i++) {
        cin >> _cluster_labels[i];
    }

    auto ret = checkConvexForCV(_grid_asses, _cluster_labels);
    cout << ret[0] << ", " << ret[1] << endl;
}

int main() {
    test2();
    return 0;
}