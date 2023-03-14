#ifndef _MEASURE_ALPHA_TRIPLES_H
#define _MEASURE_ALPHA_TRIPLES_H

/************************************
 * 
 * Author: Jiashu
 * Date: 2023-03-14
 * Description: Convexity measure of part triples: A New Convexity Measure Based on a Probabilistic Interpretation of Images
 *      part triples by alpha division
 * 
 * ***********************************/

#include "../utils/base.h"
#include "../utils/geometry.h"

static std::vector<double> grid_pos(int point){
    int x = point / grid_square;
    int y = point % grid_square;
    return std::vector<double> {1.0 * x, 1.0 * y};
}

std::vector<double> checkConvexForAlphaTriples(
    const std::vector<int> &grid_asses,
    const std::vector<int> &cluster_labels,
    double alpha = 0.5)
{
    // get position of each grid cell and calculate
    grid_square = round(sqrt(grid_asses.size()));
    std::vector<std::vector<double>> points;
    for(auto i = 0; i < grid_asses.size();i ++) {
        points.push_back(grid_pos(i));
    }
    int max_label = 0;
    for(auto i = 0; i < grid_asses.size();i++) max_label = max(max_label, cluster_labels[i] + 1);
    std::vector<std::vector<double>> label_pairs(max_label, std::vector<double>(2, 0));
    for(auto i = 0; i < grid_asses.size();i ++) {
        for(auto j = 0; j < grid_asses.size(); j++) {
            if(i == j) continue;
            if(cluster_labels[grid_asses[i]] != cluster_labels[grid_asses[j]]) continue;
            auto label = cluster_labels[grid_asses[i]];
            auto & pos1 = points[i];
            auto & pos2 = points[j];
            int alpha_grid = grid(round(pos1[0] * alpha + pos2[0] *(1 - alpha)), round(pos1[1] * alpha + pos2[1] * (1 - alpha)));
            if(cluster_labels[alpha_grid] == label) {
                label_pairs[label][0] += 1;
            }
            label_pairs[label][1] += 1;
        }
    }
    double x = 0, y = 0;
    for(auto & pair: label_pairs) {
        x += pair[0];
        y += pair[1];
    }
    return {x, y};
}

#endif