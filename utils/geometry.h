#ifndef _GEOMETRY_H
#define _GEOMETRY_H

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>

void getSurroundRect(
    double x1, double y1, double x2, double y2,
    double &xmin, double &xmax, double &ymin, double &ymax 
) {
    if(x1 < x2) {
        xmin = x1; xmax = x2;
    }
    else {
        xmin = x2; xmax = x1;
    }
    if(y1 < y2) {
        ymin = y1; ymax = y2;
    }
    else {
        ymin = y2; ymax = y1;
    }
}

bool judgeSegmentsIntersect(
    double x1, double y1, double x2, double y2,
    double x3, double y3, double x4, double y4)
{
    double xmin1, xmax1, ymin1, ymax1;
    double xmin2, xmax2, ymin2, ymax2;
    getSurroundRect(x1, y1, x2, y2, xmin1, xmax1, ymin1, ymax1);
    getSurroundRect(x3, y3, x4, y4, xmin2, xmax2, ymin2, ymax2);
	if(xmax1 < xmin2 || xmax2 < xmin1 ||ymax1 < ymin2 || ymax2 < ymin1)   return false;

	int res1 = (x1 - x2) * (y3 - y2) - (y1 - y2) * (x3 - x2);
	int res2 = (x1 - x2) * (y4 - y2) - (y1 - y2) * (x4 - x2);
	
	int res3 = (x3 - x4) * (y1 - y4) - (y3 - y4) * (x1 - x4);
	int res4 = (x3 - x4) * (y2 - y4) - (y3 - y4) * (x2 - x4);
	if(res1 * res2 <= 0 && res3 * res4 <= 0) return true;
	else return false;
}

double getAngle(const std::vector<double> & vc1, const std::vector<double> vc2) {
    // calculate Angle in two vectors [-pi, pi]
    double dot = vc1[0] * vc2[0] + vc1[1] * vc2[1];
    double l1 = sqrt(vc1[0] * vc1[0] + vc1[1] * vc1[1]);
    double l2 = sqrt(vc2[0] * vc2[0] + vc2[1] * vc2[1]);
    double angle = acos(dot / (l1 * l2));
    // std::cout << angle << " "<< l1 << " "<<l2 <<"\n";
    double skew_product = vc1[0] * vc2[1] - vc1[1] * vc2[0];
    if(skew_product >= 0) return angle;
    else return -1 * angle;
}

double getPerimeter(const std::vector<std::vector<double>> & vcs, const int type = 0) {
    double res = 0;
    if(type == 0) {
        for(auto i = 1;i < vcs.size(); i++) {
            auto & vc1 = vcs[i];
            auto & vc2 = vcs[i-1]; 
            res += sqrt((vc1[0] - vc2[0]) * (vc1[0] - vc2[0]) + (vc1[1] - vc2[1]) * (vc1[1] - vc2[1]));
        }
        auto & vc1 = vcs[0];
        auto & vc2 = vcs[vcs.size()-1]; 
        res += sqrt((vc1[0] - vc2[0]) * (vc1[0] - vc2[0]) + (vc1[1] - vc2[1]) * (vc1[1] - vc2[1]));
    }
    else {
        for(auto & vc: vcs) {
            res += sqrt(vc[0] * vc[0] + vc[1] * vc[1]);
        }
    }
    return res;
}


#endif