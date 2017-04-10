#include "utils.h"
#include "segment.h"
#include "RandomClass.h"
#include <iostream>


point emit_from_triangle(point pt1, point pt2, point pt3,RandomClass * r)
{
    bool found = 0;
    point pt;
    double R1;
    double R2;
    while (!found){
        R1 = r->randu();
        R2 = r->randu();
        if (R1+R2<1)
            found = 1;
    }
    pt.x = pt1.x+R1*(pt2.x-pt1.x)+R2*(pt3.x-pt1.x);
    pt.y = pt1.y+R1*(pt2.y-pt1.y)+R2*(pt3.y-pt1.y);
    return pt;
}


int choose(RandomClass * r, double norm_cumul[],int N_interv)
{
	//assign a frequency
	double R=r->randu();
	int i1 = 0;
	int i2 = (N_interv+1)/2;
	int i3 = N_interv+1;
	int count = 0;
	while (i2-i1>0){
	    count=count+1;
	    if (R<norm_cumul[i2-1]){
	        i3 = i2;
	        i2 = (i2+i1)/2;}
	    else{
	        i1 = i2;
	        i2 = (i2+i3)/2;
	    }
		//    cout << norm_cumul[i2-1] << ": " << i2 << endl;

	}
	//	cout << R << " " << norm_cumul[i2-1] << ": " << i2 << endl;
	return i2;
}

double det(segment seg1, segment seg2)
{
    return det(seg1.point2.x-seg1.point1.x, seg1.point2.y - seg1.point1.y,seg2.point2.x-seg2.point1.x, seg2.point2.y - seg2.point1.y);
}

point emit_from_segment(segment seg,RandomClass * r)
{
    point pt;
    double R1 = r->randu();
    pt.x = seg.point1.x+R1*(seg.point2.x-seg.point1.x);
    pt.y = seg.point1.y+R1*(seg.point2.y-seg.point1.y);
    return pt;
}
