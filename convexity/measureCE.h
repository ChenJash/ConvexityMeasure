#ifndef _MEASURE_CE_H
#define _MEASURE_CE_H


#include "../utils/base.h"
#include "../utils/convexHull.h"
#include "../utils/util.h"
#include "measureCHC.h"
#include "measureEdge.h"

//calculate the full cost, using convexity based on perimeter(where C means Circumference) of convex hull
//save: if save convexity measure pairs of every cluster
//load: if load saved convexity measure pairs of every cluster and only re-calculate clusters of "mainLabel1" and "mainLabel2"
//mainLabel1, mainLabel2: re-calculate clusters of "mainLabel1" and "mainLabel2" when load saved convexity measure pairs
//return cost[0..3], full/Similar/Compact/Convex cost
std::vector<double> checkCostForCE(
    const double Similar_cost_matrix[],
    const double Compact_cost_matrix[],
    const int grid_asses[], const int cluster_labels[],
    const int &N, const int &num, const int &square_len, const int &maxLabel,
    const double &alpha, const double &beta,
    bool save=false, bool load=false, double label_pairs[][2]=nullptr, int mainLabel1=-1, int mainLabel2=-1) {

    std::vector<double> E_pair = checkConvexForEArray(grid_asses, cluster_labels, N, num, square_len, maxLabel);
    std::vector<double> C_pair = checkConvexForCArray(grid_asses, cluster_labels, N, num, square_len, maxLabel, save, load, label_pairs, mainLabel1, mainLabel2);
    double correct_E=0, full_E=0, correct_C=0, full_C=0;
    full_E = E_pair[1];
    correct_E = E_pair[1]-E_pair[0];
    full_C = C_pair[1];
    correct_C = C_pair[1]-C_pair[0];

    double Convex_cost = (full_E-correct_E)/full_E+(full_C-correct_C)/full_C;
    double Similar_cost = 0;
    double Compact_cost = 0;
    double cost = 0;

    int *element_asses = new int[num];
    for(int i=0;i<N;i++)if(grid_asses[i]<num)element_asses[grid_asses[i]] = i;

    for(int i=0;i<num;i++){
        Similar_cost += Similar_cost_matrix[element_asses[i]*N+i];
        Compact_cost += Compact_cost_matrix[element_asses[i]*N+i];
        cost += Similar_cost_matrix[element_asses[i]*N+i] * (1-beta-alpha);
        cost += Compact_cost_matrix[element_asses[i]*N+i] * beta;
    }
    cost += Convex_cost * N * alpha;

    delete[] element_asses;

    std::vector<double> ret(4, 0);
    ret[0] = cost;
    ret[1] = Similar_cost;
    ret[2] = Compact_cost;
    ret[3] = Convex_cost*N;
    return ret;
}

#endif