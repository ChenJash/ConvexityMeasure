#ifndef _MEASURE_ACOK_H
#define _MEASURE_ACOK_H

/************************************
 * 
 * Author: Jiashu
 * Date: 2023-03-08
 * Description: Convexity measure of AocK from paper: A Measure of Non-convexity in the Plane and the Minkowski Sum
 * 
 * ***********************************/

#include "../utils/base.h"
#include "../utils/geometry.h"

typedef std::vector<double> Point, Vector;
typedef std::vector<std::vector<double>> PointList, VectorList;

double checkPolylineAocK(const VectorList & polyline) {
    // return aocK of polyline with double (-infinite, 0]
    std::vector<double> angles;
    double angle_sum = 0;
    for(auto i = 1; i < polyline.size(); i++) {
        double angle = getAngle(polyline[i-1], polyline[i]);
        angles.push_back(angle);
        angle_sum += angle;
    }
    double angle = getAngle(polyline[polyline.size()-1], polyline[0]);
    angles.push_back(angle);
    angle_sum += angle;
    if(angle_sum < 0) {
        for(auto i = 0; i < angles.size(); i++) {
            angles[i] = -1 * angles[i];
        }
        angle_sum = -1 * angle_sum;
    }

    // dp find max depress part
    std::vector<double> cur_heights;
    double left_max_height = 0, cur_height = 0, min_rot = 0;
    for(auto i = 0; i < angles.size(); i++) {
        cur_heights.push_back(cur_height);
        cur_height += angles[i];
        double cur_rot = cur_height - left_max_height;
        if(cur_rot < min_rot) min_rot = cur_rot;
        if(cur_height > left_max_height) left_max_height = cur_height;
    }
    cur_heights.push_back(cur_height);

    std::vector<int> lf(angles.size(), 0), rt(angles.size(), 0);
    double lf_min = 0;
    for(auto i = 0; i < angles.size(); i++) {
        if(cur_heights[i] < lf_min) {
            lf[i] = cur_heights[i];
            lf_min = lf[i];
        }
        else {
            lf[i] = lf_min;
        }
    }
    double rt_max = 0;
    for(auto i = angles.size() - 1; i > -1; i--) {
        double temp = cur_heights[i + 1] - angle_sum;
        if(temp > rt_max) {
            rt[i] = temp;
            rt_max = rt[i];
        }
        else {
            rt[i] = rt_max;
        }
    }
    for(auto i = 0; i < angles.size(); i++) {
        double rot = lf[i] - rt[i];
        if(rot < min_rot) {
            min_rot = rot;
        }
    }

    return min_rot;
}

std::vector<double> checkPolygonConvexByAocK(
    const PointList & polygon,
    const std::vector<PointList> & holes) 
{
    // judge polygon simple
    VectorList vectors;
    if(!judgePolygonSimple(polygon, vectors)) return std::vector<double> {0,-1};

    std::vector<VectorList> holes_vectors;
    for(auto & hole: holes) {
        VectorList hole_vectors;
        if(!judgePolygonSimple(hole, hole_vectors)) return std::vector<double> {0,-1};
        holes_vectors.push_back(hole_vectors);
    }

    // calculate rot list, get aocK by dp
    double AocK = checkPolylineAocK(vectors);
    double perimeter = getPerimeter(vectors, 1);
    for(auto hole_vectors: holes_vectors) {
        double alpha = getPerimeter(hole_vectors, 1) / perimeter;
        AocK += alpha * (checkPolylineAocK(hole_vectors) - 2 * M_PI);
    }
    // to be discussed: translate way
    return std::vector<double> {exp(AocK / M_PI), 1};
}

std::vector<double> checkConvexForAocK(
    const std::vector<int> &_grid_asses,
    const std::vector<int> &_cluster_labels)
{   
    // get boundary of clusters and check convexity of data
    std::vector<PointList> rt_boundary;
    std::vector<int> rt_boundary_labels;
    getClusterBoundary(_grid_asses, _cluster_labels, rt_boundary, rt_boundary_labels);

    // check convexity of clusters
    int cur_idx = 0;
    std::vector<int> labels;
    std::vector<double> perimeters;
    std::vector<std::vector<double>> results;
    while(cur_idx < rt_boundary.size()) {
        int cur_label = rt_boundary_labels[cur_idx];
        double max_perimeter = -1;
        int max_idx = -1;
        int fin_idx = cur_idx;
        double total_perimeter = 0;
        while(fin_idx < rt_boundary.size() && rt_boundary_labels[fin_idx] == cur_label) {
            double pm = getPerimeter(rt_boundary[fin_idx]);
            if(pm > max_perimeter) {
                max_perimeter = pm;
                max_idx = fin_idx; 
            }
            total_perimeter += pm;
            fin_idx += 1;
        }
        PointList & polygon = rt_boundary[max_idx]; 
        std::vector<PointList> holes;
        for(auto i = cur_idx; i < fin_idx; i++) {
            if(i == max_idx) continue;
            holes.push_back(rt_boundary[i]);
        }
        labels.push_back(cur_label);
        perimeters.push_back(total_perimeter);
        results.push_back(checkPolygonConvexByAocK(polygon, holes));
        cur_idx = fin_idx;
    }

    // for(auto i = 0;i < labels.size(); i++){
    //     std::cout << labels[i] << " " << perimeters[i] << " " << results[i][0] << std::endl;
    // }

    // calculate summary
    double a = 0, b = 0;
    for(auto i = 0;i < labels.size(); i++){
        a += perimeters[i] * results[i][0];
        b += perimeters[i] * results[i][1];
    }
    return std::vector<double> {a, b};
}

#endif