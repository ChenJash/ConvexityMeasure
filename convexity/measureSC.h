#ifndef _MEASURE_SC_H
#define _MEASURE_SC_H

/************************************
 * 
 * Author: Jiashu
 * Date: 2023-03-10
 * Description: Convexity measure of SC from paper: A Convexity Measure for Open and Closed Contours
 * 
 * ***********************************/

#include <algorithm>
#include "../utils/base.h"
#include "../utils/geometry.h"
#include "../utils/convexHull.h"

typedef std::vector<double> Point, Vector;
typedef std::vector<std::vector<double>> PointList, VectorList;

void mergeRepeatedLength(
    std::vector<double> & rt_lengths,
    std::vector<double> & rt_angles) 
{
    double cur_length = rt_lengths[0], cur_angle = rt_angles[0];
    int cur_idx = 0;
    for(auto i = 1; i < rt_lengths.size(); i++){
        if(fabs(rt_angles[i] - cur_angle) < 1e-5){
            cur_length += rt_lengths[i];
        }
        else{
            rt_lengths[cur_idx] = cur_length;
            rt_angles[cur_idx] = cur_angle;
            cur_idx += 1;
            cur_length = rt_lengths[i];
            cur_angle = rt_angles[i];
        }
    }
    rt_lengths[cur_idx] = cur_length;
    rt_angles[cur_idx] = cur_angle;
    cur_idx += 1;
    rt_lengths.erase(rt_lengths.begin() + cur_idx, rt_lengths.end());
    rt_angles.erase(rt_angles.begin() + cur_idx, rt_angles.end());
}

void getLengthAndAngles(
    const VectorList edge_vectors,
    std::vector<double> & rt_lengths,
    std::vector<double> & rt_angles,
    bool merge = false) 
{
    // calculate lengths and angles
    rt_lengths.clear();
    rt_angles.clear();
    double angle_sum = 0, length_sum = 0;
    double length = getVectorLength(edge_vectors[0]);
    double angle = getAngle(edge_vectors[edge_vectors.size()-1], edge_vectors[0]);
    rt_lengths.push_back(length);
    rt_angles.push_back(angle);
    length_sum += length;
    angle_sum += angle;
    for(auto i = 1; i < edge_vectors.size(); i++){
        length = getVectorLength(edge_vectors[i]);
        angle = getAngle(edge_vectors[i-1], edge_vectors[i]);
        rt_lengths.push_back(length);
        rt_angles.push_back(angle);
        length_sum += length;
        angle_sum += angle;
    }

    // normalization
    for(auto i = 0; i < rt_lengths.size(); i++){
        rt_lengths[i] = rt_lengths[i] / length_sum;
    }
    if(angle_sum < 0) {
        std::reverse(rt_lengths.begin() + 1, rt_lengths.end());
        std::reverse(rt_angles.begin() + 1, rt_angles.end());
        for(auto i = 0; i < rt_angles.size(); i++) {
            rt_angles[i] = -1 * rt_angles[i];
        }
        angle_sum = -1 * angle_sum;
    }
    double cur_height = 0;
    for(auto i = 1; i < rt_angles.size(); i++){
        rt_angles[i] = rt_angles[i] + rt_angles[i-1];
    }
    if(merge) {
        mergeRepeatedLength(rt_lengths, rt_angles);
    }
}

double checkPolylineSC(
    const PointList & polyline,
    VectorList & edge_vectors) 
{
    // get polygon convex hull
    double (*polyline_cp)[2] = new double [polyline.size()][2];
    for(auto i = 0;i < polyline.size(); i++) {
        polyline_cp[i][0] = polyline[i][0];
        polyline_cp[i][1] = polyline[i][1];
    }
    int num_nodes = getConvexHull(polyline.size(), polyline_cp);
    VectorList hull_vectors;
    Point start {polyline_cp[0][0], polyline_cp[0][1]};
    for(auto i = 1;i < num_nodes; i++){
        hull_vectors.push_back(Vector{polyline_cp[i][0]-polyline_cp[i-1][0], polyline_cp[i][1] - polyline_cp[i-1][1]});
    }
    hull_vectors.push_back(Vector{polyline_cp[0][0]-polyline_cp[num_nodes-1][0], polyline_cp[0][1] - polyline_cp[num_nodes-1][1]});
    delete [] polyline_cp;

    // get the same start point
    int start_pos = 0;
    static double epsilon = 1e-5;
    for(;start_pos < polyline.size(); start_pos++){
        if(fabs(polyline[start_pos][0] - start[0]) < epsilon && 
            fabs(polyline[start_pos][1] - start[1]) < epsilon)
            break;
    }
    if(start_pos == polyline.size()) {
        std::cout << "Convex hull calculate error." << std::endl;
        return 0;
    }
    std::rotate(edge_vectors.begin(), edge_vectors.begin() + start_pos, edge_vectors.end());


    // calulate rot angles functions and merge
    std::vector<double> ln_raw, ln_hull;
    std::vector<double> angle_raw, angle_hull;
    getLengthAndAngles(edge_vectors, ln_raw, angle_raw, true);
    getLengthAndAngles(hull_vectors, ln_hull, angle_hull);

    std::vector<double> ln_merge, angle_merge_raw, angle_merge_hull;
    double cur_pos = 0;
    double cur_id_raw = 0, cur_id_hull = 0;
    double nxt_pos_raw = ln_raw[0], nxt_pos_hull = ln_hull[0];
    while(fabs(cur_pos - 1) > epsilon) {
        if(cur_id_hull < ln_hull.size() && nxt_pos_hull < nxt_pos_raw) {
            ln_merge.push_back(nxt_pos_hull - cur_pos);
            angle_merge_raw.push_back(angle_raw[cur_id_raw]);
            angle_merge_hull.push_back(angle_hull[cur_id_hull]);
            cur_id_hull += 1;
            cur_pos = nxt_pos_hull;
            if(cur_id_hull < ln_hull.size()) nxt_pos_hull += ln_hull[cur_id_hull];
        } else if(cur_id_raw < ln_raw.size() && nxt_pos_raw < nxt_pos_hull) {
            ln_merge.push_back(nxt_pos_raw - cur_pos);
            angle_merge_raw.push_back(angle_raw[cur_id_raw]);
            angle_merge_hull.push_back(angle_hull[cur_id_hull]);
            cur_id_raw += 1;
            cur_pos = nxt_pos_raw;
            if(cur_id_raw < ln_raw.size()) nxt_pos_raw += ln_raw[cur_id_raw];
        } else if(cur_id_hull < ln_hull.size() && cur_id_raw < ln_raw.size() && nxt_pos_hull == nxt_pos_raw){
            ln_merge.push_back(nxt_pos_raw - cur_pos);
            angle_merge_raw.push_back(angle_raw[cur_id_raw]);
            angle_merge_hull.push_back(angle_hull[cur_id_hull]);
            cur_id_hull += 1;
            cur_id_raw += 1;
            cur_pos = nxt_pos_raw;
            if(cur_id_hull < ln_hull.size()) nxt_pos_hull += ln_hull[cur_id_hull];
            if(cur_id_raw < ln_raw.size()) nxt_pos_raw += ln_raw[cur_id_raw];
        }
    }

    // calulate best args \theta
    double theta = 0;
    for(auto i = 0; i < ln_merge.size(); i++) {
        theta -= (angle_merge_raw[i] - angle_merge_hull[i]) * ln_merge[i];
    }
    // calulate DH and SC
    double DH = 0;
    for(auto i = 0; i < ln_merge.size(); i++) {
        DH += pow(angle_merge_raw[i] - angle_merge_hull[i] + theta, 2) * ln_merge[i];
    }
    return 1 / (1 + DH);
}

std::vector<double> checkPolygonConvexBySC(
    const PointList & polygon,
    const std::vector<PointList> & holes) 
{   
    // check polygon simple
    VectorList vectors;
    if(!judgePolygonSimple(polygon, vectors)) return std::vector<double> {0,-1};

    std::vector<VectorList> holes_vectors;
    for(auto & hole: holes) {
        VectorList hole_vectors;
        if(!judgePolygonSimple(hole, hole_vectors)) return std::vector<double> {0,-1};
        holes_vectors.push_back(hole_vectors);
    }

    // check polyline SC and get total SC
    double punish_factor = 0.5;
    double perimeter = getPerimeter(vectors, 1);
    double SC = checkPolylineSC(polygon, vectors) * perimeter;
    double total_perimeter = perimeter;
    for(auto i = 0;i < holes.size(); i++) {
        double hole_perimeter = getPerimeter(holes_vectors[i], 1);
        SC += checkPolylineSC(holes[i], holes_vectors[i]) * hole_perimeter * punish_factor;
        total_perimeter += hole_perimeter;
    }
    SC = SC / total_perimeter;
    return std::vector<double> {SC, 1};
}

std::vector<double> checkConvexForSC(
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
        results.push_back(checkPolygonConvexBySC(polygon, holes));
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
