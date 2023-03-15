#ifndef _MEASURE_CV_H
#define _MEASURE_CV_H

#include "../utils/convexHull.h"
#include <iostream>
#include <float.h>
#include <cmath>

const double tolerance = 0.0000001; // constant for comparing two `double`s

enum Direction {
    GO_UP,
    GO_DOWN,
    GO_LEFT,
    GO_RIGHT
};

int getLabel(
const std::vector<std::vector<int>> &grid_label,
const int &x, const int &y) {
    int n = grid_label.size();
    if (x < 0 || y < 0 || x >= n || y >= n) return -1;
    return grid_label[x][y];
}

double getCV(
const int grids[][2],
const int grids_n) {
    // pre-requirements:
    // grids[] are ordered counter-clockwisely
    double (*convex_hull)[2] = new double[grids_n][2];
    for (int i = 0; i < grids_n; i++) {
        convex_hull[i][0] = grids[i][0];
        convex_hull[i][1] = grids[i][1];
    }
    int vertex_num = getConvexHull(grids_n, convex_hull);

    // L: the perimeter of convex hull
    double L = 0;
    for (int i = 0; i < vertex_num; i++) {
        int j = (i + 1) % vertex_num;
        double x1 = convex_hull[i][0];
        double y1 = convex_hull[i][1];
        double x2 = convex_hull[j][0];
        double y2 = convex_hull[j][1];
        L += std::sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    }
    // std::cout << "L: " << L << '\n';

    // 找到 grids[] 里与 convex_hull[0] 相同的顶点
    int index = 0;
    for (int i = 0; i < grids_n; i++) {
        if (std::abs(grids[i][0] - convex_hull[0][0]) < tolerance 
            && std::abs(grids[i][1] - convex_hull[0][1]) < tolerance) {
            index = i;
            break;
        }
    }

    double temp = 0;
    for (int pos = 0; pos < vertex_num; pos++) {
        int new_pos = (pos + 1) % vertex_num;
        double x1 = convex_hull[pos][0];
        double y1 = convex_hull[pos][1];
        double x2 = convex_hull[new_pos][0];
        double y2 = convex_hull[new_pos][1];

        int index2 = (index + 1) % grids_n;
        if (std::abs(grids[index2][0]-x2) < tolerance
        && std::abs(grids[index2][1]-y2) < tolerance) {
            index = index2;
            continue;
        }

        temp -= std::sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

        // debug
        // double lbi = std::sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
        // std::cout << "Lbi: " << lbi << '\n';
        
        while (true) {
            int new_index = (index + 1) % grids_n;

            double rx1 = grids[index][0];
            double ry1 = grids[index][1];
            double rx2 = grids[new_index][0];
            double ry2 = grids[new_index][1];
            temp += std::sqrt((rx1-rx2)*(rx1-rx2)+(ry1-ry2)*(ry1-ry2));

            // debug
            // double li = std::sqrt((rx1-rx2)*(rx1-rx2)+(ry1-ry2)*(ry1-ry2));
            // std::cout << "Li: " << li << '\n';

            index = new_index;

            if (std::abs(grids[index][0]-x2) < tolerance
            && std::abs(grids[index][1]-y2) < tolerance) {
                break;
            }
        }
    }

    delete [] convex_hull;
    return 1 - temp / L;
}

std::vector<double> checkConvexForCVArray(
const int grid_asses[],
const int cluster_labels[],
const int &N, const int &num, const int &square_len, const int &maxLabel) {

    std::vector<std::vector<int>> grid_label(square_len);
    for (int i = 0; i < square_len; i++) {
        grid_label[i].resize(square_len);
    }
    int *head = new int[maxLabel];
    int *last = new int[num];
    int *element_asses = new int[num];
    int *label_count = new int[maxLabel];
    for (int i = 0; i < maxLabel; i++) {
        head[i] = -1;
        label_count[i] = 0;
    }
    for (int gid = 0; gid < N; gid++) {
        int x = gid % square_len;
        int y = gid / square_len;
        if (grid_asses[gid] < num) {
            int id = grid_asses[gid];
            int lb = cluster_labels[id];
            element_asses[id] = gid;
            last[id] = head[lb]; head[lb] = id;
            grid_label[x][y] = lb;
            label_count[lb]++;
        } else {
            grid_label[x][y] = -1;
        }
    }

    double S0 = 0, S1 = 0;
    int (*nodes)[2] = new int[num*4][2];
    for (int li = 0; li < maxLabel; li++) {
        //std::cout << "label: " << li << '\n';
        int cnt = 0;
        double tmp_S0 = 0, tmp_S1 = 0;
        int px = square_len; int py = square_len;
        for (int id = head[li]; id >= 0; id = last[id]) {
            int gid = element_asses[id];
            int x = gid % square_len;
            int y = gid / square_len;
            if ((y < py) || (y == py && x < px)) {
                px = x;
                py = y;
            }
        }
        nodes[0][0] = px; nodes[0][1] = py; cnt = 1;
        int nx = px+1; int ny = py; int dir = GO_RIGHT;
        // 规定走向为逆时针
        while (true) {
            //std::cout << "nx: " << nx << ", ny: " << ny << '\n';
            if (nx == px && ny == py) break;
            switch (dir) {
                case GO_UP: {
                    int label = getLabel(grid_label, nx, ny);
                    if (label == li) {
                        nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                        nx++;
                        dir = GO_RIGHT;
                    } else {
                        label = getLabel(grid_label, nx-1, ny);
                        if (label != li) {
                            nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                            nx--;
                            dir = GO_LEFT;
                        } else {
                            ny++;
                            dir = GO_UP;
                        }
                    }
                    break;
                }
                case GO_DOWN: {
                    int label = getLabel(grid_label, nx-1, ny-1);
                    if (label == li) {
                        nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                        nx--;
                        dir = GO_LEFT;
                    } else {
                        label = getLabel(grid_label, nx, ny-1);
                        if (label != li) {
                            nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                            nx++;
                            dir = GO_RIGHT;
                        } else {
                            ny--;
                            dir = GO_DOWN;
                        }
                    }
                    break;
                }
                case GO_LEFT: {
                    int label = getLabel(grid_label, nx-1, ny);
                    if (label == li) {
                        nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                        ny++;
                        dir = GO_UP;
                    } else {
                        label = getLabel(grid_label, nx-1, ny-1);
                        if (label != li) {
                            nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                            ny--;
                            dir = GO_DOWN;
                        } else {
                            nx--;
                            dir = GO_LEFT;
                        }
                    }
                    break;
                }
                case GO_RIGHT: {
                    int label = getLabel(grid_label, nx, ny-1);
                    if (label == li) {
                        nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                        ny--;
                        dir = GO_DOWN;
                    } else {
                        label = getLabel(grid_label, nx, ny);
                        if (label != li) {
                            nodes[cnt][0] = nx; nodes[cnt][1] = ny; cnt++;
                            ny++;
                            dir = GO_UP;
                        } else {
                            nx++;
                            dir = GO_RIGHT;
                        }
                    }
                    break;
                }
            }
        }
        tmp_S1 = label_count[li];

        // debug
        // std::cout << "polygon:\n";
        // for (int i = 0; i < cnt; i++) {
        //     std::cout << nodes[i][0] << ", " << nodes[i][1] << '\n';
        // }

        tmp_S0 = getCV(nodes, cnt) * tmp_S1;
        std::cout << "tmp_S0: " << tmp_S0 << ", tmp_S1: " << tmp_S1 << '\n';

        S1 += tmp_S1;
        S0 += tmp_S0;
    }
    std::vector<double> S_pair(2, 0);
    S_pair[0] = S1 - S0;
    S_pair[1] = S1;

    delete [] head;
    delete [] last;
    delete [] element_asses;
    delete [] nodes;
    delete [] label_count;
    return S_pair;
}

std::vector<double> checkConvexForCV(
const std::vector<int> &_grid_asses,
const std::vector<int> &_cluster_labels) {
    int N = _grid_asses.size();
    int num = _cluster_labels.size();
    int square_len = ceil(sqrt(N));
    int maxLabel = 0;
    for(int i=0;i<num;i++)maxLabel = std::max(maxLabel, _cluster_labels[i]+1);

    int *grid_asses = new int[N];
    int *cluster_labels = new int[num];
    for(int i=0;i<N;i++)grid_asses[i] = _grid_asses[i];
    for(int i=0;i<num;i++)cluster_labels[i] = _cluster_labels[i];

    std::vector<double> ret = checkConvexForCVArray(grid_asses, cluster_labels, N, num, square_len, maxLabel);

    delete[] grid_asses;
    delete[] cluster_labels;
    
    return ret;
}

#endif