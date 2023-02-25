#ifndef _CONVEXHULL_H
#define _CONVEXHULL_H

#include <iostream>
#include <vector>
#include <utility>
#include <ctime>
#include <algorithm>
#include <math.h>

double Distance(const double nodes[][2], int a, int b)
{
    return std::sqrt((nodes[a][0]-nodes[b][0])*(nodes[a][0]-nodes[b][0])+(nodes[a][1]-nodes[b][1])*(nodes[a][1]-nodes[b][1]));
}

double Multiply(const double nodes[][2], int a, int b, int c)
{
    return((nodes[b][0]-nodes[a][0])*(nodes[c][1]-nodes[a][1])-(nodes[b][1]-nodes[a][1])*(nodes[c][0]-nodes[a][0]));
}

double (*global_nodes)[2];
int global_nodes_N = -1;

bool cmp(const int &x, const int &y){
    double m;
    m = Multiply(global_nodes, 0, x, y);
    if(m>0)return true;
    else if(m==0&&(Distance(global_nodes, 0,x)<Distance(global_nodes, 0, y)))
        return true;
    else return false;
}

// N: num of vertexes
// nodes[0..N-1][2]: vertex
// return int M，it's the num of vertexes of convex hull
// nodes will be modified，nodes[0..M-1][2] is the vertexes of convex hull in clockwise order
int getConvexHull(int N, double nodes[][2]){
    if(N<=2)return N;

    if(global_nodes_N!=N){
        if(global_nodes_N>=0){
            delete[] global_nodes;
            global_nodes_N=-1;
        }
        global_nodes = new double[N+1][2];
        global_nodes_N = N;
    }
    
    double px, py;
    int p;
    py=-1;
    for(int i=0;i<N;i++)
    {
        if(py==-1||nodes[i][1]<py)
        {
            px=nodes[i][0];
            py=nodes[i][1];
            p=i;
        }
        else if(nodes[i][1]==py&&nodes[i][0]<px)
        {
            px=nodes[i][0];
            py=nodes[i][1];
            p=i;
        }
    }
    double tmp=nodes[0][0]; nodes[0][0]=nodes[p][0]; nodes[p][0]=tmp;
    tmp=nodes[0][1]; nodes[0][1]=nodes[p][1]; nodes[p][1]=tmp;

    for(int i=0;i<N;i++){
        global_nodes[i][0] = nodes[i][0];
        global_nodes[i][1] = nodes[i][1];
    }
    global_nodes[N][0] = global_nodes[0][0]; global_nodes[N][1] = global_nodes[0][1];

    int* order = new int[N+1];
    int* order2 = new int[N+1];
    for(int i=0;i<N+1;i++)order[i]=i;
    std::sort(order+1, order+N, cmp);

    order2[0] = order[0];
    order2[1] = order[1];
    order2[2] = order[2];

    int top=2;
    for(int ii=3;ii<=N;ii++)
    {
        while(top>=1&&Multiply(global_nodes, order2[top-1], order2[top], order[ii])<=0)
            top--;
        order2[top+1] = order[ii];
        top++;
    }

    for(int i=0;i<top;i++){
        nodes[i][0] = global_nodes[order2[i]][0];
        nodes[i][1] = global_nodes[order2[i]][1];
    }

    delete[] order;
    delete[] order2;
    return top;
}

// get the area of convex hull
double getSofPoly(int N, double nodes[][2]){
    if(N<=2)return 0;
    double area = 0;
    for(int i=0;i<N-1;i++){
        double triArea = (nodes[i][0]*nodes[i+1][1] - nodes[i+1][0]*nodes[i][1])/2;
        area += triArea;
    }
    double fn = (nodes[N-1][0]*nodes[0][1] - nodes[0][0]*nodes[N-1][1])/2;
    return abs(area+fn);
}

// get the perimeter of convex hull
double getCofPoly(int N, double nodes[][2]){
    if(N<=1)return 0;
    double C = 0;
    for(int i=0;i<N-1;i++){
        double edge = Distance(nodes, i, i+1);
        C += edge;
    }
    double edge = Distance(nodes, N-1, 0);
    return abs(C+edge);
}

#endif