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

int intrsct_test(segment segmentab, segment segment12) //tests whether two ORIENTED segments intersect
{ //important: segmentab is a trajectory segment; segment12 is a boundary (any type);
  //return +1 if the collision is "direct" (determinant of vectors>0), -1 if indirect, 0 if no collision
	double dtr = det(segmentab, segment12);

	double t, s; // represent normalized coordinates on the two segments
	if (dtr!=0)
	{
		t = ((segmentab.point1.y - segment12.point1.y)*(segment12.point2.x - segment12.point1.x)
			 -(segmentab.point1.x - segment12.point1.x)*(segment12.point2.y - segment12.point1.y))/dtr;
		s = -((segment12.point1.y - segmentab.point1.y)*(segmentab.point2.x - segmentab.point1.x)
			 -(segment12.point1.x - segmentab.point1.x)*(segmentab.point2.y - segmentab.point1.y))/dtr;
		if (0<=t && t<=1 && 0<=s && s<=1)
		{
			if (dtr>0)
			{	return 1;}else{return -1;}

		}
		else {
			return 0;
		}

	}
	else {
		return 0;
	}

}

point intrsct_pt(segment segmentab, segment segment12) //calculates the intersection point between two ORIENTED segments intersect
{ //important: needs to be used along with intrsct_pt to check if the intersection exists
	double dtr = det(segmentab, segment12);
	point pt;
	double t, s; // represent normalized coordinates on the two segments

	if (dtr!=0)
	{
		t = ((segmentab.point1.y - segment12.point1.y)*(segment12.point2.x - segment12.point1.x)
			 -(segmentab.point1.x - segment12.point1.x)*(segment12.point2.y - segment12.point1.y))/dtr;
		s = -((segment12.point1.y - segmentab.point1.y)*(segmentab.point2.x - segmentab.point1.x)
			  -(segment12.point1.x - segmentab.point1.x)*(segmentab.point2.y - segmentab.point1.y))/dtr;

		pt.x = segmentab.point1.x+(segmentab.point2.x-segmentab.point1.x)*t;
		pt.y = segmentab.point1.y+(segmentab.point2.y-segmentab.point1.y)*t;
	}
	return pt;

}

double dist_to_intrsct(segment segmentab,point pt){ // calculates the distance between starting point of segment and a point (to be used with previous functions above)
	return sqrt((pt.x-segmentab.point1.x)*(pt.x-segmentab.point1.x)+(pt.y-segmentab.point1.y)*(pt.y-segmentab.point1.y));
}
