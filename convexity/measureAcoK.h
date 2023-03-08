#ifndef _MEASURE_ACOK_H
#define _MEASURE_ACOK_H

/************************************
 * 
 * Author: Jiashu
 * Date: 2023-03-08
 * Description: Convexity measure of AocK from paper: A Measure of Non-convexity in the Plane and the Minkowski Sum
 * 
 * ***********************************/
// // for test without base
// #include<iostream>
// #include<vector>
// #include<unordered_set>
// #include<unordered_map>
// #include<cmath>
// template <typename T>
// inline T max(const T&a, const T& b) {
//     return a > b ? a 
//     : b;
// }
// template <typename T>
// inline T min(const T&a, const T& b) {
//     return a < b ? a : b;
// }
// // end for test

#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include "../utils/base.h"
#include "../utils/geometry.h"

typedef std::vector<double> Point, Vector;
typedef std::vector<std::vector<double>> PointList, VectorList;

struct pair_hash
{
    template<class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const
    {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

bool judgePolygonSimple(
    const PointList & polygon,
    VectorList & rt_vectors,
    bool judge_selfinter = false)
{
    // judge non-coincide, non-reverse
    rt_vectors.clear();
    Vector before {0, 0};
    for (auto i = 1; i < polygon.size(); i++) {
        Vector temp {polygon[i][0] - polygon[i-1][0], polygon[i][1] - polygon[i-1][1]};
        if(temp[0] == 0 && temp[1] == 0) {
            std::cout << "[AocK]: judgePolygonSimple - polygon coincides." << std::endl;
            return false;
        }
        if((temp[0] * before[1] == temp[1] * before[0]) && temp[0] * before[0] < 0) {
            std::cout << "[AocK]: judgePolygonSimple - polygon reverses." << std::endl;
            return false;
        }
        rt_vectors.push_back(temp);
        before = temp;
    }
    rt_vectors.push_back(Vector {polygon[0][0] - polygon[polygon.size()-1][0], polygon[0][1] - polygon[polygon.size()-1][1]});

    // judge simple: has no self-intersections
    if(judge_selfinter) {
        for (auto i = 2; i < polygon.size(); i++) {
            for (auto j = 1; j < i - 1; j++) {
                if(judgeSegmentsIntersect(polygon[i][0], polygon[i][1], polygon[i-1][0], polygon[i-1][1],
                    polygon[j][0], polygon[j][1], polygon[j-1][0], polygon[j-1][1])) {
                    // std::cout<< polygon[i][0]<<" " <<polygon[i][1]<<" "<<polygon[i-1][0]<<" " <<polygon[i-1][0]<<std::endl;
                    // std::cout<< polygon[j][0]<<" " <<polygon[j][1]<<" "<<polygon[j-1][0]<<" " <<polygon[j-1][0]<<std::endl;
                    std::cout << "[AocK]: judgePolygonSimple - polygon intersects." << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

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
        // to be discussed: hole punish
        AocK += alpha * (checkPolylineAocK(vectors) - 2 * M_PI);
    }
    // to be discussed: translate way
    return std::vector<double> {exp(AocK / M_PI), 1};
}

static int grid_square = 1;
static int grid(int x, int y) {
    return x * grid_square + y;
}
static int point_square = 1;
static int point(int x, int y) {
    return x * point_square + y;
}
static std::vector<double> pos(int point){
    // std::cout << point <<" "<<point_square << " ";
    int x = point / point_square;
    int y = point % point_square;
    // std::cout << x <<" "<<y << " "<< 1.0 *x << " "<< 1.0 *y << std::endl;
    return std::vector<double> {1.0 * x, 1.0 * y};
}

void getClusterBoundary(
    const std::vector<int> & grid_asses,
    const std::vector<int> & cluster_labels,
    std::vector<PointList> & rt_boundary,
    std::vector<int> & rt_boundary_labels)
{
    // get cluster boundary segs
    int cur_id = 0;
    std::unordered_map<int, int> label2ids;
    std::vector<std::unordered_set<std::pair<int, int>, pair_hash>> edges;
    std::vector<std::unordered_map<int, std::vector<int>>> matrix;
    int square_len = ceil(sqrt(grid_asses.size()));
    grid_square = square_len;
    point_square = square_len + 1;

    for(auto i = 0;i < grid_asses.size(); i++) {
        int gid = grid_asses[i];
        int label = cluster_labels[gid];
        if(label2ids.find(label) == label2ids.end()){
            label2ids[label] = cur_id;
            edges.push_back(std::unordered_set<std::pair<int, int>, pair_hash>());
            matrix.push_back(std::unordered_map<int, std::vector<int>>());
            cur_id += 1;
        }
        int lid = label2ids[label];
        auto & ledges = edges[lid];
        auto & lmatrix = matrix[lid];
        
        int x = i / square_len;
        int y = i % square_len;
        int edge_start, edge_end;
        if(x == 0 || cluster_labels[grid_asses[grid(x-1, y)]] != label) {
            edge_start = point(x, y);
            edge_end = point(x, y+1);
            ledges.insert(std::pair<int, int>(edge_start, edge_end));
            lmatrix[edge_start].push_back(edge_end);
            lmatrix[edge_end].push_back(edge_start);
        }
        if(y == 0 || cluster_labels[grid_asses[grid(x, y-1)]] != label) {
            edge_start = point(x, y);
            edge_end = point(x+1, y);
            ledges.insert(std::pair<int, int>(edge_start, edge_end));
            lmatrix[edge_start].push_back(edge_end);
            lmatrix[edge_end].push_back(edge_start);
        }
        if(x == square_len-1 || cluster_labels[grid_asses[grid(x+1, y)]] != label) {
            edge_start = point(x+1, y);
            edge_end = point(x+1, y+1);
            ledges.insert(std::pair<int, int>(edge_start, edge_end));
            lmatrix[edge_start].push_back(edge_end);
            lmatrix[edge_end].push_back(edge_start);
        }
        if(y == square_len-1 || cluster_labels[grid_asses[grid(x, y+1)]] != label) {
            edge_start = point(x, y+1);
            edge_end = point(x+1, y+1);
            ledges.insert(std::pair<int, int>(edge_start, edge_end));
            lmatrix[edge_start].push_back(edge_end);
            lmatrix[edge_end].push_back(edge_start);
        }
    }

    // connect boundary
    rt_boundary.clear();
    rt_boundary_labels.clear();
    for(auto it = label2ids.begin(); it != label2ids.end(); it++) {
        int label = it->first;
        int lid = it->second;
        auto & ledges = edges[lid];
        auto & lmatrix = matrix[lid];
        while(!ledges.empty()) {
            auto edge_head = *ledges.begin();
            ledges.erase(ledges.begin());
            PointList points;
            points.push_back(pos(edge_head.first));
            int start = edge_head.first;
            int now = edge_head.second;
            int bf = start;
            while(now != start) {
                auto & link = lmatrix[now];
                int nxt = 0;
                if(link[0] == bf) {
                    nxt = link[1];
                }
                else nxt = link[0];
                ledges.erase(std::pair<int, int>(min(now, nxt), max(now, nxt)));
                points.push_back(pos(now));
                bf = now;
                now = nxt;
            }
            rt_boundary.push_back(points);
            rt_boundary_labels.push_back(label);
        }
    }
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