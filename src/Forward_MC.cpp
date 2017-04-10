/*
Copyright (c) 2016, Jean-Philippe M. Péraud
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS AND CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/



// a set of classes and functions:
// classes:
// - segment: stores coordinates of two points and calculate the length
// - reflective: inherited from segment, also has a member spec indicating the specularity of a reflective boundary
// - prescribed: inherited from segment, also has a member Teq (equilibrium temperature) and T (actual temperature)
// - periodic: ////////////////////////, also has a member vctr for indicating the translation to be applied to a particle
// - quadrilater: a set of four segments with area and center.
// - reflective_bdrs: reads and stores all reflective boundaries
// - prescribed_bdrs: reads and stores all prescribed temperature boundaries
// - periodic_bdrs: reads and stores all periodic boundaries
// - materials: stores the material properties
// - sources

// functions:
// -det: calculated determinant between two vectors
// - calc_area: calculates area of triangle defined by three points
// - intrsct_test: determines whether two segments intersect
// - intrsct_pt: determines intersection point between to segments
// - dist_to_inrsct: determines distance between starting point of segment and another specified point
// - inside_quad: tests whether a point is inside a quadrilater
// - overlap_quad: test whether a segment overlaps with a quadrilater

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdio.h>
#include <cstring>
#include <string>
#include <ctime>
#include <math.h>
#include "randomClass.h"

using namespace std;
double PI = 3.14159265;

struct point{
    double x;
    double y;
};

double det(double ax, double ay, double bx, double by)
{
	if (ax*by-ay*bx>1) {
	//	cout << "in det: too large: " << ax << " " << ay << " " << bx << " " << by << endl;
	}
    return ax*by-ay*bx;
}

double calc_area(point pt1, point pt2, point pt3)
{
    return det(pt2.x-pt1.x, pt2.y-pt1.y, pt3.x-pt1.x, pt3.y-pt1.y)/2;
};

double dist_pts(point pt1, point pt2)
{
	return sqrt((pt2.x-pt1.x)*(pt2.x-pt1.x)+(pt2.y-pt1.y)*(pt2.y-pt1.y));
}

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

class segment {
public:
    segment(){};
    segment(double x1, double y1, double x2, double y2){point1.x = x1; point1.y = y1; point2.x = x2; point2.y = y2; length = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));}
    point point1;
    point point2;
    double length;
	point normal;
    void assign(double x1, double y1, double x2, double y2){point1.x = x1; point1.y = y1; point2.x = x2; point2.y = y2;}
	void initialize(); // calculates normal and length if not already done
	void show(){cout << "x1 : " << point1.x << " y1 : " << point1.y << " x2 : " << point2.x << " y2 : " << point2.y << endl;}

};

void segment::initialize()
{
	double x1, x2, y1, y2;
	x1 = point1.x; y1 = point1.y; x2 = point2.x; y2 = point2.y;
	length = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    if (length!=0)
    {
        normal.x = -(y2-y1)/length;
        normal.y = (x2-x1)/length;
    }
    else{
        normal.x = 0;
        normal.y = 0;
    }
}

double det(segment seg1, segment seg2)
{
    return det(seg1.point2.x-seg1.point1.x, seg1.point2.y - seg1.point1.y,seg2.point2.x-seg2.point1.x, seg2.point2.y - seg2.point1.y);
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////



class reflective: public segment {
public:
	reflective(){};
    reflective(double x1, double y1, double x2, double y2, double sp);
    double spec; //specularity coefficient (1 = totally specular)
};

reflective::reflective(double x1, double y1, double x2, double y2, double sp){
    point1.x = x1; point1.y = y1; point2.x = x2; point2.y = y2; spec = sp;
    length = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    if (length!=0)
    {
        normal.x = -(y2-y1)/length;
        normal.y = (x2-x1)/length;
    }
    else{
        normal.x = 0;
        normal.y = 0;
    }
}


class prescribed: public segment {
public:
	prescribed(){};
    prescribed(double x1, double y1, double x2, double y2, double T, double Teq);
    double temp; //this is the DEVIATIONAL temperature
    double temp_eq;
};

prescribed::prescribed(double x1, double y1, double x2, double y2, double T, double Teq){
    point1.x = x1; point1.y = y1; point2.x = x2; point2.y = y2; temp = T; temp_eq = Teq;
    length = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    if (length!=0)
    {
        normal.x = -(y2-y1)/length;
        normal.y = (x2-x1)/length;
    }
    else{
        normal.x = 0;
        normal.y = 0;
    }
}

class periodic: public segment{
public:
    periodic(){};
    periodic(double x1, double y1, double x2, double y2, double shiftx, double shifty);
    point vctr; // this stores the two translation coordinates to be applied to a particle when it hits this wall
};

periodic::periodic(double x1, double y1, double x2, double y2, double shiftx, double shifty){
    point1.x = x1; point1.y = y1; point2.x = x2; point2.y = y2; vctr.x = shiftx; vctr.y = shifty;
    length = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    if (length!=0)
    {
        normal.x = -(y2-y1)/length;
        normal.y = (x2-x1)/length;
    }
    else{
        normal.x = 0;
        normal.y = 0;
    }
}

point emit_from_segment(segment seg,RandomClass * r)
{
    point pt;
    double R1 = r->randu();
    pt.x = seg.point1.x+R1*(seg.point2.x-seg.point1.x);
    pt.y = seg.point1.y+R1*(seg.point2.y-seg.point1.y);
    return pt;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////






class quadrilater{ //class that describes a 2D region, calculates and stores its area
	// segments/points must be organized such that it rotates in the direct sense
public:
    quadrilater(){};
    quadrilater(segment segment1, segment segment2, segment segment3, segment segment4);
    quadrilater(point pt1, point pt2, point pt3, point pt4);
    segment segs[4];
    double area;
    point center;
};

quadrilater::quadrilater(segment segment1, segment segment2, segment segment3, segment segment4)
{
    // assign points to constitutive segments
    segs[0].point1.x = segment1.point1.x;
    segs[0].point1.y = segment1.point1.y;
    segs[0].point2.x = segment1.point2.x;
    segs[0].point2.y = segment1.point2.y;

    segs[1].point1.x = segment2.point1.x;
    segs[1].point1.y = segment2.point1.y;
    segs[1].point2.x = segment2.point2.x;
    segs[1].point2.y = segment2.point2.y;

    segs[2].point1.x = segment3.point1.x;
    segs[2].point1.y = segment3.point1.y;
    segs[2].point2.x = segment3.point2.x;
    segs[2].point2.y = segment3.point2.y;

    segs[3].point1.x = segment4.point1.x;
    segs[3].point1.y = segment4.point1.y;
    segs[3].point2.x = segment4.point2.x;
    segs[3].point2.y = segment4.point2.y;

    // calculate area
    area = calc_area(segs[0].point1,segs[1].point1,segs[1].point2)
         +calc_area(segs[0].point1,segs[2].point1,segs[2].point2);

    // center of the shape
    center.x = (segment1.point1.x + segment2.point1.x+segment3.point1.x + segment4.point1.x)/4;
    center.y = (segment1.point1.y + segment2.point1.y+segment3.point1.y + segment4.point1.y)/4;

}

quadrilater::quadrilater(point pt1, point pt2, point pt3, point pt4)
{
    // assign points to constitutive segments
    segs[0].point1 = pt1;
    segs[0].point2 = pt2;

    segs[1].point1 = pt2;
    segs[1].point2 = pt3;

    segs[2].point1 = pt3;
    segs[2].point2 = pt4;

    segs[3].point1 = pt4;
    segs[3].point2 = pt1;

    // calculate area
    area = calc_area(segs[0].point1,segs[1].point1,segs[1].point2)
         +calc_area(segs[0].point1,segs[2].point1,segs[2].point2);

    // center of the shape
    center.x = (pt1.x + pt2.x + pt3.x + pt4.x)/4;
    center.y = (pt1.y + pt2.y + pt3.y + pt4.y)/4;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////








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

bool inside_quad(point pt, quadrilater quad) // checks if point pt is inside quad. IMPORTANT: only works if quad is convex
{
	// just test by checking that following determinants are all positive
	bool test = 1;
	for (int i = 0; i<4; i++){
		if (det(quad.segs[i].point1.x-pt.x, quad.segs[i].point1.y-pt.y, quad.segs[i].point2.x-pt.x, quad.segs[i].point2.y-pt.y)<=0){test = 0;}
	}
	return test;
}

bool overlap_quad(segment seg, quadrilater quad) // test whether a segment and a quadrilater overlap
{
	bool test = 0;
	if (inside_quad(seg.point1, quad) || inside_quad(seg.point2,quad)) // if one the two points is inside, then there is obviously some overlap
	{
		test = 1;
	}
	else {
		for (int i = 0; i<4; i++){
			if (intrsct_test(seg, quad.segs[i])!=0){test = 1;} // if the segment intersects with one of the quad side, then there is overlap
		}
	}
	return test;

}

double overlap_length(segment seg, quadrilater quad) // calculates the length of segment overlap
{// IMPORTANT: if this function is used although there is no overlap, it can return a non zero value. The function overlap_quad
// MUST be used before using this function.
	double length = 0;
	bool found1 = 0;
	bool found2 = 0;
	int i;
	point pt1;
	point pt2;
	// needs to distinguish between several cases

	if (inside_quad(seg.point1, quad) && inside_quad(seg.point2, quad)){ // both points inside

		length = sqrt((seg.point2.x-seg.point1.x)*(seg.point2.x-seg.point1.x)+(seg.point2.y-seg.point1.y)*(seg.point2.y-seg.point1.y));
	}
	if (inside_quad(seg.point1, quad) && !inside_quad(seg.point2, quad)){ // only starting point inside

		i = -1;
		while (!found1){
			i++;
			if (intrsct_test(seg, quad.segs[i])==1)
				found1 = 1;
		}
		pt2 = intrsct_pt(seg, quad.segs[i]);
		length = sqrt((pt2.x-seg.point1.x)*(pt2.x-seg.point1.x)+(pt2.y-seg.point1.y)*(pt2.y-seg.point1.y));
	}
	if (!inside_quad(seg.point1, quad) && inside_quad(seg.point2, quad)){ // only ending point inside

		i = -1;
		while (!found1){
			i++;
			if (intrsct_test(seg, quad.segs[i])==-1)
				found1 = 1;
		}

		pt1 = intrsct_pt(seg, quad.segs[i]);
		length = sqrt((seg.point2.x-pt1.x)*(seg.point2.x-pt1.x)+(seg.point2.y-pt1.y)*(seg.point2.y-pt1.y));
	}
	if (!inside_quad(seg.point1, quad) && !inside_quad(seg.point2, quad)){ // no point inside => two intersection points
		i = -1;
		while (!found1 && i<4){
			i++;
			if (intrsct_test(seg, quad.segs[i])==-1) // looking for entry point
				found1 = 1;
		}
	/*	cout <<seg.point1.x << " " << seg.point1.y << " " << seg.point2.x << " "<< seg.point2.y <<endl;
		quad.segs[0].show();
		quad.segs[1].show();
		quad.segs[2].show();
		quad.segs[3].show();
		cout << "index 1 " << i << endl;*/
		if (found1) {
		    pt1 = intrsct_pt(seg, quad.segs[i]);
		}

//		cout << "pt1 " << pt1.x << " " << pt1.y << endl;
		found2 = 0;
		i = -1;
		while (!found2 && i<4){
			i++;
			if (intrsct_test(seg, quad.segs[i])==1) // looking for exit point
				found2 = 1;
		}
//		cout << "index " << i << endl;
		if (found1) {
		    pt2 = intrsct_pt(seg, quad.segs[i]);
		}

//		cout << "pt2 " << pt2.x << " " << pt2.y << endl;
		if (found1 && found2) {
		    length = sqrt((pt2.x-pt1.x)*(pt2.x-pt1.x)+(pt2.y-pt1.y)*(pt2.y-pt1.y));
		}
		else {
			length = 0;
		}

	}

	return length;
}

point emit_from_quadrilater(quadrilater quad, RandomClass * r)
{
    // subdivide into 2 triangles and calculate their respective aboslute areas
    double A1 = abs(calc_area(quad.segs[0].point1,quad.segs[1].point1,quad.segs[1].point2));
    double A2 = abs(calc_area(quad.segs[0].point1,quad.segs[2].point1,quad.segs[2].point2));
    double threshold = A1/(A1+A2);
    double R1 = r->randu();
	point pt;
    if (R1<threshold)
    {
        pt = emit_from_triangle(quad.segs[0].point1,quad.segs[1].point1,quad.segs[1].point2, r);

    }
    else
    {
        pt = emit_from_triangle(quad.segs[0].point1,quad.segs[2].point1,quad.segs[2].point2, r);
    }
	return pt;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////






class reflective_bdrs { // class that handles and stores the reflective boundaries
public:
	reflective_bdrs(const char* filename);
	~reflective_bdrs();
	int N; // TOTAL number of reflective boundaries
	reflective * ref_handle; // pointer towards an array of reflective boundaries
	void show(); //displays all reflective boundaries
};

reflective_bdrs::reflective_bdrs(const char* filename)
{
	ifstream file;
	file.open(filename);
	if (file.fail())
	{
		cout << "no input file for reflective boundaries" << endl;
		N=0;
		ref_handle = NULL;
	}
	else {
		file >> N;
		ref_handle = new reflective[N];
		for (int i=0; i<N; i++){
			file >> ref_handle[i].point1.x;
			file >> ref_handle[i].point1.y;
			file >> ref_handle[i].point2.x;
			file >> ref_handle[i].point2.y;
			file >> ref_handle[i].spec;
			ref_handle[i].initialize();
		}
	}
}

void reflective_bdrs::show()
{
	for (int i = 0; i<N ; i++){
		cout << ref_handle[i].point1.x << " " << ref_handle[i].point1.y << " "
		    << ref_handle[i].point2.x << " " << ref_handle[i].point2.x << " "
		    << ref_handle[i].spec << " " << ref_handle[i].length << endl;
	}
}

reflective_bdrs::~reflective_bdrs()
{
	delete ref_handle;
}

class prescribed_bdrs { // class that handles and stores the prescribed temperature boundaries
public:
	prescribed_bdrs(const char* filename);
	~prescribed_bdrs();
	int N; // TOTAL number of reflective boundaries
	prescribed * presc_handle; // pointer towards an array of prescribed boundaries
	void show(); //displays all reflective boundaries
};

prescribed_bdrs::prescribed_bdrs(const char* filename)
{
	ifstream file;
	file.open(filename);
	if (file.fail())
	{
		cout << "no input file for prescribed temperature boundaries" << endl;
		N=0;
		presc_handle = NULL;
	}
	else {
		file >> N;
		presc_handle = new prescribed[N];
		for (int i=0; i<N; i++){
			file >> presc_handle[i].point1.x;
			file >> presc_handle[i].point1.y;
			file >> presc_handle[i].point2.x;
			file >> presc_handle[i].point2.y;
			file >> presc_handle[i].temp;
			file >> presc_handle[i].temp_eq;
			presc_handle[i].initialize();
		}
	}
}

void prescribed_bdrs::show()
{
	for (int i = 0; i<N ; i++){
		cout << presc_handle[i].point1.x << " " << presc_handle[i].point1.y << " "
		<< presc_handle[i].point2.x << " " << presc_handle[i].point2.x << " "
		<< presc_handle[i].temp << " " << presc_handle[i].temp_eq << endl;
	}
}

prescribed_bdrs::~prescribed_bdrs()
{
	delete presc_handle;
}

class periodic_bdrs { // class that handles and stores the prescribed temperature boundaries
public:
	periodic_bdrs(const char* filename);
	~periodic_bdrs();
	int N; // TOTAL number of reflective boundaries
	periodic * per_handle; // pointer towards an array of reflective boundaries
	void show(); //displays all reflective boundaries
};

periodic_bdrs::periodic_bdrs(const char* filename)
{
	ifstream file;
	file.open(filename);
	if (file.fail())
	{
		cout << "no input file for periodic boundaries" << endl;
		N=0;
		per_handle = NULL;
	}
	else {
		file >> N;
		per_handle = new periodic[N];
		for (int i=0; i<N; i++){
			file >> per_handle[i].point1.x;
			file >> per_handle[i].point1.y;
			file >> per_handle[i].point2.x;
			file >> per_handle[i].point2.y;
			file >> per_handle[i].vctr.x;
			file >> per_handle[i].vctr.y;
			per_handle[i].initialize();
		}
	}
}

void periodic_bdrs::show()
{
	for (int i = 0; i<N ; i++){
		cout << per_handle[i].point1.x << " " << per_handle[i].point1.y << " "
		<< per_handle[i].point2.x << " " << per_handle[i].point2.x << " "
		<< per_handle[i].vctr.x << " " << per_handle[i].vctr.y << endl;
	}
}

periodic_bdrs::~periodic_bdrs()
{
	delete per_handle;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////






class materials
{
public:
        materials(const char* filename, double Teq); // initialize with file of material properties
        ~materials();
        // including frequencies, densities of states, velocities, relaxation times, Delta frequencies, polarizations.
        // IMPORTANT NOTE: the factor 2 of the TA modes must already be included in the corresponding densities of state.
	int Nm; //total number of modes
	double C; //heat capacity
	double * F; //frequencies
	double * SD; // density of states
	double * VG; //group velocities
	double * tau; // relaxation times
	double * Dom; // delta of frequencies
	int * pol;

	double hbar; //reduced planck constant
	double boltz; //Boltzmann constant
	double * de_dT; // derivative of Bose Einstein

	// useful distributions
	double * C_om;
	double * Ctau;
	double * CV;
	// cumulative distributions (for source terms)
	double * cumul_C;
	double * cumul_Ctau;
	double * cumul_CV;
	// normalized cumulative distributions
	double * N_cumul_C;
	double * N_cumul_Ctau;
	double * N_cumul_CV;

	double kth; // thermal conductivity

	void show_all(); // displays all distributions and data

};

materials::materials(const char* filename, double Teq)
{ // filename must obey very strict data organization, as follows:
// 1st entry = total number of "modes" (or "frequency bins")
// then each line = frequency, density of states, velocity, Delta frequency, relaxation time, polarization
	// Teq is the equilibrium temperature, necessary for calculating the temperature dependent properties
	hbar = 1.054517e-34;
	boltz = 1.38065e-23;
	ifstream file;
	file.open(filename);
	if (file.fail())
	{
		cout << "no input file for material data" << endl;
		Nm=0;
		F = NULL;
		SD = NULL;
		VG = NULL;
		tau = NULL;
		Dom = NULL;
		pol = NULL;
		de_dT = NULL;
		C_om = NULL;
		Ctau = NULL;
		CV = NULL;
		cumul_C = NULL;
		cumul_Ctau = NULL;
		cumul_CV = NULL;
		N_cumul_C = NULL;
		N_cumul_Ctau = NULL;
		N_cumul_CV = NULL;
	}
	else {
		file >> Nm;
		F = new double[Nm];
		SD = new double[Nm];
		VG = new double[Nm];
		tau = new double[Nm];
		Dom = new double[Nm];
		pol = new int[Nm];
		de_dT = new double[Nm];
		C_om = new double[Nm];
		Ctau = new double[Nm];
		CV = new double[Nm];
		cumul_C = new double[Nm];
		cumul_Ctau = new double[Nm];
		cumul_CV = new double[Nm];
		N_cumul_C = new double[Nm];
		N_cumul_Ctau = new double[Nm];
		N_cumul_CV = new double[Nm];
		// read material parameters
		for (int i=0; i<Nm; i++){
			file >> F[i];
			file >> SD[i];
			file >> VG[i];
			file >> Dom[i];
			file >> tau[i];
			file >> pol[i];
			de_dT[i] = (hbar*F[i]/Teq)*(hbar*F[i]/Teq)/boltz*exp(hbar*F[i]/boltz/Teq)/(exp(hbar*F[i]/boltz/Teq)-1)/(exp(hbar*F[i]/boltz/Teq)-1);
			C_om[i] = SD[i]*de_dT[i]*Dom[i];
			Ctau[i] = SD[i]*de_dT[i]*Dom[i]/tau[i];
			CV[i] = SD[i]*de_dT[i]*Dom[i]*VG[i];
		}

		// calculate cumulative distributions
		cumul_C[0] = C_om[0];
		cumul_Ctau[0] = Ctau[0];
		cumul_CV[0] = CV[0];
		kth = CV[0]*VG[0]*tau[0]/3;
		for (int i=1; i<Nm; i++){
			cumul_C[i] = cumul_C[i-1] + C_om[i];
			cumul_Ctau[i] = cumul_Ctau[i-1] + Ctau[i];
			cumul_CV[i] = cumul_CV[i-1] + CV[i];
			kth = kth + CV[i]*VG[i]*tau[i]/3;
		}

		// normalize cumulative distributions
		for (int i = 0; i<Nm; i++)
		{
			N_cumul_C[i] = cumul_C[i]/cumul_C[Nm-1];
			N_cumul_Ctau[i] = cumul_Ctau[i]/cumul_Ctau[Nm-1];
			N_cumul_CV[i] = cumul_CV[i]/cumul_CV[Nm-1];
		}
		C = cumul_C[Nm-1];
	}
}

materials::~materials()
{
	delete[] F;
	delete[] SD;
	delete[] VG;
	delete[] tau;
	delete[] Dom;
	delete[] pol;
	delete[] de_dT;
	delete[] C_om;
	delete[] Ctau;
	delete[] CV;
	delete[] cumul_C;
	delete[] cumul_Ctau;
	delete[] cumul_CV;
	delete[] N_cumul_C;
	delete[] N_cumul_Ctau;
	delete[] N_cumul_CV;
}

void materials::show_all()
{
        cout << "number of frequency cells: " << Nm << endl;
	cout << "Boltzmann constant: " << boltz << endl;
	cout << "reduced Planck constant: " << hbar << endl;
	cout << "RAW DATA: " << endl;
	for (int i = 0; i<Nm; i++){
		cout << F[i] << " " << SD[i] << " " << VG[i] << " " << Dom[i] << " " << tau[i] << " " << pol[i] << endl;
	}
	cout << "DISTRIBUTIONS: " << endl;
	for (int i = 0; i<Nm; i++){
		cout << i << " " <<C_om[i] << " " << Ctau[i] << " " << CV[i] << endl;
	}
	cout << "CUMULATIVE DISTRIBUTIONS: " << endl;
	for (int i = 0; i<Nm; i++){
		cout << i << " " <<cumul_C[i] << " " << cumul_Ctau[i] << " " << cumul_CV[i] << endl;
	}
	cout << "NORMALIZED CUMULATIVE DISTRIBUTIONS:" << endl;
	for (int i = 0; i<Nm; i++){
		cout << i << " " << N_cumul_C[i] << " " << N_cumul_Ctau[i] << " " << N_cumul_CV[i] << endl;
	}
	cout << "thermal conductivity : " << kth << endl;

}



 class particle
 {
 public:

	 particle(){COLL_MAX = 10;}
	 particle(int Collm){COLL_MAX = Collm;}
         int sig; // sign of particle, 1 or -1
	 int mode_index;
         point pt0;
         point pt1;
         double tau; //current relaxation time
         int pol; // current polarization
	 double V; //norm of velocity
         point Vp0; // velocity (2D vector)
         point Vp1;
         double t; // particle's internal time.
	 double Dt; // records the time to next scattering event
         double weight; //particle's weight
         int counter; // counts the number of times the particle has moved
	 int COLL_MAX;
         bool alive; // indicates whether the particle is still active or terminated
         segment seg; // segment used for processing trajectories and collisions with boundaries
         int collision_type; // 0= sca1 = reflection; 2 = prescribed (absorption); 3 = periodic --- useful for function "finish_move".
	 int collision_index; // to record which boundary was hit
         void initiate_move(RandomClass* r); // starts the move without looking at boundaries, and assigns seg
         void move(reflective_bdrs * ref, prescribed_bdrs * presc, periodic_bdrs * per); // uses the information in seg to locate the target point after considering boundaries
        // UPDATES seg (needed for detectors)

         void finish_move(materials * mat, RandomClass * r,reflective_bdrs * ref, prescribed_bdrs * presc, periodic_bdrs * per);
	 void show_all();
         bool P_inside_quad(quadrilater quad);
 };


void particle::show_all()
{
	cout << "x0: " << pt0.x << "  y0: " << pt0.y << "  index: " << mode_index << "  sig: " << sig <<
	 "  pol: " << pol << "  V: " << V << "  Vp0x: " << Vp0.x << "  Vp0y: " << Vp0.y << "  t: " << t << "  Dt: " << Dt <<
	 "  x1: " << pt1.x << "  y1: " << pt1.y << " last collision type:  " << collision_type <<endl;
}

void particle::initiate_move(RandomClass * r)
{
	Dt = -log(1-r->randu())*tau;
	pt1.x = pt0.x + Dt*Vp0.x;
	pt1.y = pt0.y + Dt*Vp0.y;
	seg.point1.x = pt0.x;
	seg.point1.y = pt0.y;
	seg.point2.x = pt1.x;
	seg.point2.y = pt1.y;
	seg.length = Dt*sqrt(Vp0.x*Vp0.x+Vp0.y*Vp0.y); // DO NOT USE V BECAUSE 2D
}

void particle::move(reflective_bdrs * ref, prescribed_bdrs * presc, periodic_bdrs * per)
{
	point test_pt;
	double pdist;
	collision_type = 0;
	for (int i = 0; i<ref->N; i++) // try and find the closest boundary that is being hit, if any
	{
		if (intrsct_test(seg, ref->ref_handle[i])==1) {
			test_pt = intrsct_pt(seg, ref->ref_handle[i]);
			pdist = dist_pts(pt0,test_pt);
			if (seg.length>pdist) {
				pt1 = test_pt;
				seg.point2 = pt1;
				seg.length = pdist;
				collision_type = 1;
				collision_index = i;
				Dt = pdist/sqrt(Vp0.x*Vp0.x+Vp0.y*Vp0.y); // present version does not handle velocity 0; //!! in 2D cannot use V, must use sqrt(Vx^2+Vy^2)
			}


		}
	}
	for (int i = 0; i<presc->N; i++)
	{
		if (intrsct_test(seg, presc->presc_handle[i])==1) {
			test_pt = intrsct_pt(seg, presc->presc_handle[i]);
			pdist = dist_pts(pt0,test_pt);
			if (seg.length>pdist) {
				pt1 = test_pt;
				seg.point2 = pt1;
				seg.length = pdist;
				collision_type = 2;
				collision_index = i;
				Dt = pdist/sqrt(Vp0.x*Vp0.x+Vp0.y*Vp0.y);
			}
		}
	}
	for (int i = 0; i<per->N; i++)
	{
		if (intrsct_test(seg, per->per_handle[i])==1) {
			test_pt = intrsct_pt(seg, per->per_handle[i]);
			pdist = dist_pts(pt0,test_pt);
			if (seg.length>pdist) {
				pt1 = test_pt;
				seg.point2 = pt1;
				seg.length = pdist;
				collision_type = 3;
				collision_index = i;
				Dt = pdist/sqrt(Vp0.x*Vp0.x+Vp0.y*Vp0.y);
			}
		}
	}
}

void particle::finish_move(materials * mat, RandomClass * r, reflective_bdrs * ref, prescribed_bdrs * presc, periodic_bdrs * per)
{
    double R, phi, nx, ny, Vx, Vy;
    t = t + Dt; // new time. important for transient cases
	switch(collision_type){
		case 0  : // scattering (no boundary involved)
			// update position
			pt0 = pt1;
			// find new mode index
			mode_index = choose(r, mat->N_cumul_Ctau, mat->Nm);
			// update velocity
			V = mat->VG[mode_index];
            R = 2*r->randu()-1;
            phi = 2*PI*r->randu();
            Vp0.x = V*R;
            Vp0.y = V*sqrt(1-R*R)*cos(phi);
			// update other properties
			tau = mat->tau[mode_index]; //current relaxation time
	        pol = mat->pol[mode_index]; // current polarization
			counter++; //increment collision counter
            break;
		case 1  : // collision with reflective boundary
			// update position
            pt0 = pt1;
			// mode index is unchanged

			// update velocity
			// determine if reflection is specular or diffuse
			R = r->randu();
			if (R>ref->ref_handle[collision_index].spec) // diffuse reflection
            {
			    nx = ref->ref_handle[collision_index].normal.x; ny = ref->ref_handle[collision_index].normal.y;
                R = r->randu(); phi = 2*PI*r->randu(); Vx = sqrt(R)*V; Vy = sqrt(1-R)*V*cos(phi);
	    	    Vp0.x = nx*Vx-ny*Vy;
		        Vp0.y = ny*Vx+nx*Vy;
				counter++; //counts as actual collision
            }
            else
            { // specular reflection
                nx = ref->ref_handle[collision_index].normal.x; ny = ref->ref_handle[collision_index].normal.y;
                Vx = -Vp0.x*nx-Vp0.y*ny;
                Vy = -Vp0.x*ny+Vp0.y*nx;
                Vp0.x = nx*Vx-ny*Vy;
		        Vp0.y = ny*Vx+nx*Vy;
            }

			// other properties unchanged
            break;
		case 2  : // collision with prescribed boundary
			// simply terminate the particle
			alive = false;
		//	show_all();
            break;
		case 3  : // collision with periodic boundary
			// update position
			pt0.x = pt1.x + per->per_handle[collision_index].vctr.x;
			pt0.y = pt1.y + per->per_handle[collision_index].vctr.y;

			// all the rest is unchanged (does not count as collision, therefore counter unchanged too
			break;
	}

	if (counter>COLL_MAX){alive = false;}
}

bool particle::P_inside_quad(quadrilater quad)
{
	return inside_quad(pt0, quad);
}

class volumetric: public quadrilater
{
public:
	volumetric(){};
	volumetric(point pt1, point pt2, point pt3, point pt4, double T, double Teq);
	double temp;
	double temp_eq;

};

volumetric::volumetric(point pt1, point pt2, point pt3, point pt4, double T, double Teq)
{
    // assign points to constitutive segments
    segs[0].point1 = pt1;
    segs[0].point2 = pt2;

    segs[1].point1 = pt2;
    segs[1].point2 = pt3;

    segs[2].point1 = pt3;
    segs[2].point2 = pt4;

    segs[3].point1 = pt4;
    segs[3].point2 = pt1;

    // calculate area
    area = calc_area(segs[0].point1,segs[1].point1,segs[1].point2)
	+calc_area(segs[0].point1,segs[2].point1,segs[2].point2);

    // center of the shape
    center.x = (pt1.x + pt2.x + pt3.x + pt4.x)/4;
    center.y = (pt1.y + pt2.y + pt3.y + pt4.y)/4;

	temp = T;
	temp_eq = Teq;

}

class initial: public quadrilater // for defining initial conditions in the transient case
{
public:
	initial(){};
	initial(point pt1, point pt2, point pt3, point pt4, double T, double Teq);
	double temp;
	double temp_eq;

};

initial::initial(point pt1, point pt2, point pt3, point pt4, double T, double Teq)
{
    // assign points to constitutive segments
    segs[0].point1 = pt1;
    segs[0].point2 = pt2;

    segs[1].point1 = pt2;
    segs[1].point2 = pt3;

    segs[2].point1 = pt3;
    segs[2].point2 = pt4;

    segs[3].point1 = pt4;
    segs[3].point2 = pt1;

    // calculate area
    area = calc_area(segs[0].point1,segs[1].point1,segs[1].point2)
	+calc_area(segs[0].point1,segs[2].point1,segs[2].point2);

    // center of the shape
    center.x = (pt1.x + pt2.x + pt3.x + pt4.x)/4;
    center.y = (pt1.y + pt2.y + pt3.y + pt4.y)/4;

    temp = T;
    temp_eq = Teq;

}

class body_force: public quadrilater
{
public:
	body_force(){};
	body_force(point pt1, point pt2, point pt3, point pt4, point grad);
	point vctr;

};

body_force::body_force(point pt1, point pt2, point pt3, point pt4, point grad)
{
    // assign points to constitutive segments
    segs[0].point1 = pt1;
    segs[0].point2 = pt2;

    segs[1].point1 = pt2;
    segs[1].point2 = pt3;

    segs[2].point1 = pt3;
    segs[2].point2 = pt4;

    segs[3].point1 = pt4;
    segs[3].point2 = pt1;

    // calculate area
    area = calc_area(segs[0].point1,segs[1].point1,segs[1].point2)
	+calc_area(segs[0].point1,segs[2].point1,segs[2].point2);

    // center of the shape
    center.x = (pt1.x + pt2.x + pt3.x + pt4.x)/4;
    center.y = (pt1.y + pt2.y + pt3.y + pt4.y)/4;

    vctr = grad;

}

class detector_H: public quadrilater
{
public:
	detector_H(){estimates = NULL;};
	~detector_H();
	detector_H(point pt1, point pt2, point pt3, point pt4, point grad);
	point vctr;
	double estimate;
	double * estimates;
};

detector_H::detector_H(point pt1, point pt2, point pt3, point pt4, point grad) // grad is the "imposed" temperature gradient
{
    // assign points to constitutive segments
    segs[0].point1 = pt1;
    segs[0].point2 = pt2;

    segs[1].point1 = pt2;
    segs[1].point2 = pt3;

    segs[2].point1 = pt3;
    segs[2].point2 = pt4;

    segs[3].point1 = pt4;
    segs[3].point2 = pt1;


    vctr.x = grad.x/sqrt(grad.x*grad.x+grad.y*grad.y);
    vctr.y = grad.y/sqrt(grad.x*grad.x+grad.y*grad.y);

    // calculate area
    area = calc_area(segs[0].point1,segs[1].point1,segs[1].point2)
	+calc_area(segs[0].point1,segs[2].point1,segs[2].point2);

    // center of the shape
    center.x = (pt1.x + pt2.x + pt3.x + pt4.x)/4;
    center.y = (pt1.y + pt2.y + pt3.y + pt4.y)/4;

    estimate = 0;
}

detector_H::~detector_H()
{
	if (estimates!=NULL) {
		delete[] estimates;
	}
}

class detector_T: public quadrilater
{
public:
	detector_T(){estimates = NULL;};
	~detector_T();
	detector_T(point pt1, point pt2, point pt3, point pt4);
    double estimate;
	double * estimates; // reserved for the transient case
};

detector_T::detector_T(point pt1, point pt2, point pt3, point pt4)
{
    // assign points to constitutive segments
    segs[0].point1 = pt1;
    segs[0].point2 = pt2;

    segs[1].point1 = pt2;
    segs[1].point2 = pt3;

    segs[2].point1 = pt3;
    segs[2].point2 = pt4;

    segs[3].point1 = pt4;
    segs[3].point2 = pt1;

    // calculate area
    area = calc_area(segs[0].point1,segs[1].point1,segs[1].point2)
	+calc_area(segs[0].point1,segs[2].point1,segs[2].point2);

    // center of the shape
    center.x = (pt1.x + pt2.x + pt3.x + pt4.x)/4;
    center.y = (pt1.y + pt2.y + pt3.y + pt4.y)/4;

    estimate = 0;

}

detector_T::~detector_T()
{
	if (estimates!=NULL) {
		delete[] estimates;
	}
}

class detector_array_H { // class that handles and stores the heat flux detectors
public:
	detector_array_H(const char* filename, const char * filename_time);
	~detector_array_H();
	void measure(particle * part); // updates all heat flux detectors from particle trajectory
	int N; // TOTAL number of heat flux detectors
	int Nt;
	string type; // TRANSIENT or STEADY
	detector_H * h_handle; // pointer towards an array of heat flux detectors
	double * msr_times;
	int msr_index;
	void show(); //displays all heat flux detectors
	void show_results(); // displays all heat flux results
	void write(const char* filename); // writes all heat flux results in a file
};

detector_array_H::detector_array_H(const char* filename, const char * filename_time)
{
	ifstream file, filetime;
	file.open(filename);
	filetime.open(filename_time);
	if (filetime.fail()) {
		Nt = 0;
		type = "STEADY";
		msr_times = NULL;
	}
	else{
		type = "TRANSIENT";
		filetime >> Nt;
		msr_times = new double[Nt];
		msr_index = -1;
		for (int i = 0; i<Nt; i++)
		{
			filetime >> msr_times[i];
		}
	}
	double vx, vy;
	if (file.fail())
	{
		cout << "no input file for heat flux detectors" << endl;
		N=0;
		h_handle = NULL;
	}
	else {
		file >> N;
		h_handle = new detector_H[N];
		for (int i=0; i<N; i++){
			file >> h_handle[i].segs[0].point1.x;
			file >> h_handle[i].segs[0].point1.y;
			file >> h_handle[i].segs[1].point1.x;
			file >> h_handle[i].segs[1].point1.y;
			file >> h_handle[i].segs[2].point1.x;
			file >> h_handle[i].segs[2].point1.y;
			file >> h_handle[i].segs[3].point1.x;
			file >> h_handle[i].segs[3].point1.y;
			file >> vx;
			file >> vy;
			h_handle[i].vctr.x = vx/sqrt(vx*vx+vy*vy);
			h_handle[i].vctr.y = vy/sqrt(vx*vx+vy*vy);

			h_handle[i].segs[0].point2.x = h_handle[i].segs[1].point1.x;
			h_handle[i].segs[0].point2.y = h_handle[i].segs[1].point1.y;
			h_handle[i].segs[1].point2.x = h_handle[i].segs[2].point1.x;
			h_handle[i].segs[1].point2.y = h_handle[i].segs[2].point1.y;
			h_handle[i].segs[2].point2.x = h_handle[i].segs[3].point1.x;
			h_handle[i].segs[2].point2.y = h_handle[i].segs[3].point1.y;
			h_handle[i].segs[3].point2.x = h_handle[i].segs[0].point1.x;
			h_handle[i].segs[3].point2.y = h_handle[i].segs[0].point1.y;

			h_handle[i].area = calc_area(h_handle[i].segs[0].point1,h_handle[i].segs[1].point1,h_handle[i].segs[1].point2)
			+calc_area(h_handle[i].segs[0].point1,h_handle[i].segs[2].point1,h_handle[i].segs[2].point2);

			h_handle[i].estimate = 0;
			if (strcmp(type.c_str(),"TRANSIENT")==0){
				h_handle[i].estimates = new double[Nt];
				for (int j = 0; j<Nt; j++){
			        h_handle[i].estimates[j] = 0;
				}
			}
			else {
				h_handle[i].estimate = 0;
			}
		}
	}
}

void detector_array_H::show()
{
	for (int i = 0; i<N ; i++){
		cout << h_handle[i].segs[0].point1.x << " " << h_handle[i].segs[0].point1.y << " "
		<< h_handle[i].segs[0].point2.x << " " << h_handle[i].segs[0].point2.y << " "
		<< h_handle[i].segs[1].point1.x << " " << h_handle[i].segs[1].point1.y << " "
		<< h_handle[i].segs[1].point2.x << " " << h_handle[i].segs[1].point2.y << " "
		<< h_handle[i].segs[2].point1.x << " " << h_handle[i].segs[2].point1.y << " "
		<< h_handle[i].segs[2].point2.x << " " << h_handle[i].segs[2].point2.y << " "
		<< h_handle[i].segs[3].point1.x << " " << h_handle[i].segs[3].point1.y << " "
		<< h_handle[i].segs[3].point2.x << " " << h_handle[i].segs[3].point2.y << " "
		<< h_handle[i].vctr.x << " " << h_handle[i].vctr.y  << endl;

	}
}

void detector_array_H::show_results()
{
	for (int i = 0; i<N ; i++){
		cout << h_handle[i].estimate << endl;

	}
}


detector_array_H::~detector_array_H()
{
	delete[] h_handle;
	if (msr_times!=NULL) {
		delete[] msr_times;
	}
}

void detector_array_H::measure(particle * part)
{
	double cntrbt = 0;

	if (strcmp(type.c_str(),"TRANSIENT")==0) {
		// first, we need to spot all measurement times on the path between pt0 and pt1
		int ilm=msr_index+1;
		int inm;
		double tpositions;
		point pt;

			//        cout << Vx1 << " " << Vy1 << " " << Vz1 << endl;
/*
			if (msr_times[ilm] < part->t + part->Dt){
				inm=ilm;
				while (inm-1 < Nt && msr_times[inm+1] < part->t + part->Dt ){
					inm=inm+1;
				}
				for (int im=ilm; im<inm+1; im++){
                    // calculate the positions in the phase space at the applicable measurement times
					tpositions=msr_times[im];
					pt.x = part->pt0.x + part->Vp0.x*(tpositions-part->t);
					pt.y = part->pt0.y + part->Vp0.y*(tpositions-part->t);
                    // check whether the particle would be in any of the detectors at these moments
                    for (int i = 0; i<N; i++)
                    {
					    if (tpositions > part->t && inside_quad(pt, t_handle[i]) && im<Nt){

                            t_handle[i].estimates[im] = t_handle[i].estimates[im] + part->weight*part->sig/t_handle[i].area/mat->C;//exp(-2*dist_R2/R0/R0-beta*Zpositions)/N;

					    }

				    }
				}
				msr_index=inm;
			}
			if (part->t + part->Dt > msr_times[Nt-1])
            {
                msr_index = -1; // only for transient cases: reset the msr_index
                part->alive = 0; // "kill" the particle
            }
*/
        for (int j = 0; j<Nt; j++){
            if (msr_times[j] >=part->t && msr_times[j] < part->t+part->Dt){
                tpositions=msr_times[j];
                pt.x = part->pt0.x + part->Vp0.x*(tpositions-part->t);
                pt.y = part->pt0.y + part->Vp0.y*(tpositions-part->t);
                // check whether the particle would be in any of the detectors at these moments
                for (int i = 0; i<N; i++)
                {
                    if (inside_quad(pt, h_handle[i])){
                        h_handle[i].estimates[j] = h_handle[i].estimates[j] + part->weight*part->sig/h_handle[i].area*(h_handle[i].vctr.x*part->Vp0.x+h_handle[i].vctr.y*part->Vp0.y);
				    }

                }
            }
        }
        if (part->t + part->Dt > msr_times[Nt-1])
            {
                msr_index = -1; // only for transient cases: reset the msr_index
                part->alive = 0; // "kill" the particle
            }
	}
	else{
        point nrmlzd_seg; // point structure to store the coordinates of the normalized vector colinear with segment
	    nrmlzd_seg.x = (part->seg.point2.x - part->seg.point1.x)/
	    sqrt((part->seg.point2.x - part->seg.point1.x)*(part->seg.point2.x - part->seg.point1.x)+(part->seg.point2.y - part->seg.point1.y)*(part->seg.point2.y - part->seg.point1.y));
	    nrmlzd_seg.y = (part->seg.point2.y - part->seg.point1.y)/
	    sqrt((part->seg.point2.x - part->seg.point1.x)*(part->seg.point2.x - part->seg.point1.x)+(part->seg.point2.y - part->seg.point1.y)*(part->seg.point2.y - part->seg.point1.y));

	    for (int i = 0; i<N; i++){ // go over all H detectors and analyze interaction with particle segment
		    if (overlap_quad(part->seg, h_handle[i])) {
		        cntrbt = overlap_length(part->seg, h_handle[i]); //this just gives a length
		        cntrbt = cntrbt*(nrmlzd_seg.x*h_handle[i].vctr.x + nrmlzd_seg.y*h_handle[i].vctr.y); // projects to the desired component
		        h_handle[i].estimate = h_handle[i].estimate + part->sig*part->weight*cntrbt/h_handle[i].area;
		    }
	    }
	}
}

void detector_array_H::write(const char* filename)
{
    ofstream file;
	file.open(filename);
	for (int i = 0; i<N ; i++){
        if (strcmp(type.c_str(),"TRANSIENT")==0){
            for (int j = 0; j<Nt; j++)
            {
                file << h_handle[i].estimates[j];
                file << " " ;
            }
            file << endl;
        }
        else{
		    file << h_handle[i].estimate << endl;
        }
	}
}


class detector_array_T { // class that handles and stores the temperature detectors
public:
	detector_array_T(const char* filename, const char * filename_time);
	~detector_array_T();
	void measure(particle * part, materials * mat); // updates all temperature detectors from particle trajectory
	int N; // TOTAL number of temperature detectors
	int Nt; // number of times for measurement
	string type;  //TRANSIENT or STEADY
	detector_T * t_handle; // pointer towards an array of heat flux detectors
	double * msr_times;
        int msr_index; // measurement index (for transient cases)
	void show(); //displays all temperature detectors
	void show_results();
	void write(const char * filename); // writes all temperature results in a file
};

detector_array_T::detector_array_T(const char* filename,const char* filename_time)
{
	ifstream file, filetime;
	file.open(filename);
	filetime.open(filename_time);
	if (filetime.fail()) {
		Nt = 0;
		type = "STEADY";
		msr_times = NULL;
	}
	else{
		type = "TRANSIENT";
		filetime >> Nt;
		msr_times = new double[Nt];
		msr_index = -1;
		for (int i = 0; i<Nt; i++)
		{
			filetime >> msr_times[i];
		}
	}
	double vx, vy;
	if (file.fail())
	{
		cout << "no input file for temperature detectors" << endl;
		N=0;
		t_handle = NULL;
	}
	else {
		file >> N;
		t_handle = new detector_T[N];
		for (int i=0; i<N; i++){
			file >> t_handle[i].segs[0].point1.x;
			file >> t_handle[i].segs[0].point1.y;
			file >> t_handle[i].segs[1].point1.x;
			file >> t_handle[i].segs[1].point1.y;
			file >> t_handle[i].segs[2].point1.x;
			file >> t_handle[i].segs[2].point1.y;
			file >> t_handle[i].segs[3].point1.x;
			file >> t_handle[i].segs[3].point1.y;

			t_handle[i].segs[0].point2.x = t_handle[i].segs[1].point1.x;
			t_handle[i].segs[0].point2.y = t_handle[i].segs[1].point1.y;
			t_handle[i].segs[1].point2.x = t_handle[i].segs[2].point1.x;
			t_handle[i].segs[1].point2.y = t_handle[i].segs[2].point1.y;
			t_handle[i].segs[2].point2.x = t_handle[i].segs[3].point1.x;
			t_handle[i].segs[2].point2.y = t_handle[i].segs[3].point1.y;
			t_handle[i].segs[3].point2.x = t_handle[i].segs[0].point1.x;
			t_handle[i].segs[3].point2.y = t_handle[i].segs[0].point1.y;

			t_handle[i].area = calc_area(t_handle[i].segs[0].point1,t_handle[i].segs[1].point1,t_handle[i].segs[1].point2)
			+calc_area(t_handle[i].segs[0].point1,t_handle[i].segs[2].point1,t_handle[i].segs[2].point2);
                        if (strcmp(type.c_str(),"TRANSIENT")==0){
				t_handle[i].estimates = new double[Nt];
				for (int j = 0; j<Nt; j++){
			        t_handle[i].estimates[j] = 0;
				}
			}
			else {
				t_handle[i].estimate = 0;
			}

		}
	}
}

void detector_array_T::measure(particle * part, materials * mat)
{
	double cntrbt = 0;

        if (strcmp(type.c_str(),"TRANSIENT")==0) {
		// first, we need to spot all measurement times on the path between pt0 and pt1
		int ilm=msr_index+1;
		int inm;
		double tpositions;
		point pt;
			//        cout << Vx1 << " " << Vy1 << " " << Vz1 << endl;
/*
			if (msr_times[ilm] < part->t + part->Dt){
				inm=ilm;
				while (inm-1 < Nt && msr_times[inm+1] < part->t + part->Dt ){
					inm=inm+1;
				}
				for (int im=ilm; im<inm+1; im++){
                    // calculate the positions in the phase space at the applicable measurement times
					tpositions=msr_times[im];
					pt.x = part->pt0.x + part->Vp0.x*(tpositions-part->t);
					pt.y = part->pt0.y + part->Vp0.y*(tpositions-part->t);
                    // check whether the particle would be in any of the detectors at these moments
                    for (int i = 0; i<N; i++)
                    {
					    if (tpositions > part->t && inside_quad(pt, t_handle[i]) && im<Nt){

                            t_handle[i].estimates[im] = t_handle[i].estimates[im] + part->weight*part->sig/t_handle[i].area/mat->C;//exp(-2*dist_R2/R0/R0-beta*Zpositions)/N;

					    }

				    }
				}
				msr_index=inm;
			}
			if (part->t + part->Dt > msr_times[Nt-1])
                        {
                                msr_index = -1; // only for transient cases: reset the msr_index
                                part->alive = 0; // "kill" the particle
                        }
*/
                for (int j = 0; j<Nt; j++){
                        if (msr_times[j] >=part->t && msr_times[j] < part->t+part->Dt){
                        tpositions=msr_times[j];
                        pt.x = part->pt0.x + part->Vp0.x*(tpositions-part->t);
                        pt.y = part->pt0.y + part->Vp0.y*(tpositions-part->t);
                        // check whether the particle would be in any of the detectors at these moments
                        for (int i = 0; i<N; i++)
                        {
                                if (inside_quad(pt, t_handle[i])){

                                        t_handle[i].estimates[j] = t_handle[i].estimates[j] + part->weight*part->sig/t_handle[i].area/mat->C;

			        }

                        }
                }
        }
        if (part->t + part->Dt > msr_times[Nt-1])
                {
                        msr_index = -1; // only for transient cases: reset the msr_index
                        part->alive = 0; // "kill" the particle
                 }
	}
	else{
	        for (int i = 0; i<N; i++){ // go over all T detectors and analyze interaction with particle segment
		        if (overlap_quad(part->seg, t_handle[i])) {
		                cntrbt = overlap_length(part->seg, t_handle[i]); //this just gives a length
		                t_handle[i].estimate = t_handle[i].estimate + part->sig*part->weight*
		    	        cntrbt/t_handle[i].area/mat->C/sqrt(part->Vp0.x*part->Vp0.x+part->Vp0.y*part->Vp0.y); // IMPORTANT: SINCE IT IS 2D, DO NOT USE part->V for norm of velocity
		        }
	        }
	}
}

void detector_array_T::show()
{
	for (int i = 0; i<N ; i++){
		cout << t_handle[i].segs[0].point1.x << " " << t_handle[i].segs[0].point1.y << " "
		<< t_handle[i].segs[0].point2.x << " " << t_handle[i].segs[0].point2.x << " " << t_handle[i].area
		 << endl;
	}
}

void detector_array_T::write(const char * filename)
{
	ofstream file;
	file.open(filename);
	for (int i = 0; i<N ; i++){
        if (strcmp(type.c_str(),"TRANSIENT")==0){
            for (int j = 0; j<Nt; j++)
            {
                file << t_handle[i].estimates[j];
                file << " " ;
            }
            file << endl;
        }
        else{
		    file << t_handle[i].estimate << endl;
        }
	}
}

void detector_array_T::show_results()
{
	for (int i = 0; i<N ; i++){
		cout << t_handle[i].estimate << endl;

	}
}

detector_array_T::~detector_array_T()
{
	delete[] t_handle;
	if (msr_times!=NULL) {
		delete[] msr_times;
	}
}

class sources
{
public:
	sources(materials * mat, prescribed_bdrs * presc, const char * vol, const char * bod, const char * init, int NN); // vol and bod are respectively the name files of volumetric sources and body force sources
	~sources();
	string type; //TRANSIENT or STEADY
	double total_energy;
	double t_max;  // maximum time (for TRANSIENT only)
	double * energies; // individual weights of each source
	double * cumul_energies; // same, but cumulative
	double * N_cumul; // same, but normalized
	int * source_type; // for determining from which source to emit 1 = prescribed temperature wall; 2 = volumetric heating; 3 = body force induced by temperature gradient
	prescribed_bdrs * ptr_to_presc;
	volumetric * sprd_src_array; //pointer to array of volumetric sources
	body_force * bd_frc_array; //pointer to array of body force types of sources
	initial * initial_condition;

	int Np; //number of prescribed sources
	int Nv; // number of volumetric sources
	int Nb; // number of body_force sources
	int Ni; // number of sources associated with initial condition
	int Ntot; //total number of separate sources
	int Npart; // total number of particles to be simulated
	materials * ptr_to_mat;

	void emit(particle * part, RandomClass * r); // returns a particle sampled from the sources
	void display_type(){cout << type << endl;} // to check if this is transient or steady simulation
	void display_vol();
	void display_bod();
	void display_init();
};

sources::~sources(){
	delete[] energies;
	delete[] cumul_energies;
	delete[] N_cumul;
	delete[] source_type;
	delete[] sprd_src_array;
	delete[] bd_frc_array;
	delete[] initial_condition;
}

sources::sources(materials * mat, prescribed_bdrs * presc, const char * vol, const char * bod, const char * init, int NN) // constructor
{

    Npart = NN;
	//determine the type of simulation (transient of steady)
	const char * filename_t = "times.txt";
	int Ntimes;
	ptr_to_mat = mat;
	point pt1, pt2, pt3, pt4;
	ifstream file;
	file.open(filename_t);
	if (file.fail())
	{
		cout << "no input file for times: this is a steady state calculation" << endl;
		type = "STEADY";
		Ni = 0;
	}
	else {
		type = "TRANSIENT";
		cout << "found input file for times. Estimates will be provided for times: " << endl;
		file >> Ntimes ;
		for (int i=0; i<Ntimes; i++) {
			file >> t_max;
			cout << t_max << " s" <<endl;
		}
		cout << "maximum simulation time: " << t_max << " s" << endl;
	}
	const char * local_type = type.c_str();
        file.close();
	// pointer to the prescribed sources
	ptr_to_presc = presc;
	Np = presc->N;
	// read and define volumetric
	file.open(vol);
	if (file.fail()){
		Nv = 0;
		sprd_src_array = NULL;
	}
	else
        {
	    file >> Nv; //THE FIRST NUMBER SHOULD BE THE NUMBER OF SOURCES
	    sprd_src_array = new volumetric[Nv];
	    for (int i = 0; i<Nv; i++){
		        file >> sprd_src_array[i].segs[0].point1.x;
			file >> sprd_src_array[i].segs[0].point1.y;

			file >> sprd_src_array[i].segs[1].point1.x;
			file >> sprd_src_array[i].segs[1].point1.y;

			file >> sprd_src_array[i].segs[2].point1.x;
			file >> sprd_src_array[i].segs[2].point1.y;

			file >> sprd_src_array[i].segs[3].point1.x;
			file >> sprd_src_array[i].segs[3].point1.y;

			file >> sprd_src_array[i].temp;
			file >> sprd_src_array[i].temp_eq;

			sprd_src_array[i].segs[0].point2.x = sprd_src_array[i].segs[1].point1.x;
			sprd_src_array[i].segs[0].point2.y = sprd_src_array[i].segs[1].point1.y;

			sprd_src_array[i].segs[1].point2.x = sprd_src_array[i].segs[2].point1.x;
			sprd_src_array[i].segs[1].point2.y = sprd_src_array[i].segs[2].point1.y;

			sprd_src_array[i].segs[2].point2.x = sprd_src_array[i].segs[3].point1.x;
			sprd_src_array[i].segs[2].point2.y = sprd_src_array[i].segs[3].point1.y;

			sprd_src_array[i].segs[3].point2.x = sprd_src_array[i].segs[0].point1.x;
			sprd_src_array[i].segs[3].point2.y = sprd_src_array[i].segs[0].point1.y;

			// calculate area
			sprd_src_array[i].area = calc_area(sprd_src_array[i].segs[0].point1,sprd_src_array[i].segs[1].point1,sprd_src_array[i].segs[1].point2)
			+calc_area(sprd_src_array[i].segs[0].point1,sprd_src_array[i].segs[2].point1,sprd_src_array[i].segs[2].point2);

			// center of the shape
			sprd_src_array[i].center.x = (sprd_src_array[i].segs[0].point1.x + sprd_src_array[i].segs[1].point1.x
						+ sprd_src_array[i].segs[2].point1.x + sprd_src_array[i].segs[3].point1.x)/4;
			sprd_src_array[i].center.y = (sprd_src_array[i].segs[0].point1.y + sprd_src_array[i].segs[1].point1.y
										  + sprd_src_array[i].segs[2].point1.y + sprd_src_array[i].segs[3].point1.y)/4;
	    }
	}
	file.close();
	//read and define body forces
	file.open(bod);
	if (file.fail()){
		Nb = 0;
		bd_frc_array = NULL;
	}
	else
        {
	    file >> Nb; //THE FIRST NUMBER SHOULD BE THE NUMBER OF SOURCES
	    bd_frc_array = new body_force[Nb];
	    for (int i = 0; i<Nb; i++){
		    file >> bd_frc_array[i].segs[0].point1.x;
			file >> bd_frc_array[i].segs[0].point1.y;

			file >> bd_frc_array[i].segs[1].point1.x;
			file >> bd_frc_array[i].segs[1].point1.y;

			file >> bd_frc_array[i].segs[2].point1.x;
			file >> bd_frc_array[i].segs[2].point1.y;

			file >> bd_frc_array[i].segs[3].point1.x;
			file >> bd_frc_array[i].segs[3].point1.y;

			file >> bd_frc_array[i].vctr.x;
			file >> bd_frc_array[i].vctr.y;

			bd_frc_array[i].segs[0].point2.x = bd_frc_array[i].segs[1].point1.x;
			bd_frc_array[i].segs[0].point2.y = bd_frc_array[i].segs[1].point1.y;

			bd_frc_array[i].segs[1].point2.x = bd_frc_array[i].segs[2].point1.x;
			bd_frc_array[i].segs[1].point2.y = bd_frc_array[i].segs[2].point1.y;

			bd_frc_array[i].segs[2].point2.x = bd_frc_array[i].segs[3].point1.x;
			bd_frc_array[i].segs[2].point2.y = bd_frc_array[i].segs[3].point1.y;

			bd_frc_array[i].segs[3].point2.x = bd_frc_array[i].segs[0].point1.x;
			bd_frc_array[i].segs[3].point2.y = bd_frc_array[i].segs[0].point1.y;

			// calculate area
			bd_frc_array[i].area = calc_area(bd_frc_array[i].segs[0].point1,bd_frc_array[i].segs[1].point1,bd_frc_array[i].segs[1].point2)
			+calc_area(bd_frc_array[i].segs[0].point1,bd_frc_array[i].segs[2].point1,bd_frc_array[i].segs[2].point2);

			// center of the shape
			bd_frc_array[i].center.x = (bd_frc_array[i].segs[0].point1.x + bd_frc_array[i].segs[1].point1.x
										  + bd_frc_array[i].segs[2].point1.x + bd_frc_array[i].segs[3].point1.x)/4;
			bd_frc_array[i].center.y = (bd_frc_array[i].segs[0].point1.y + bd_frc_array[i].segs[1].point1.y
										  + bd_frc_array[i].segs[2].point1.y + bd_frc_array[i].segs[3].point1.y)/4;
	    }
	}
	file.close();
    //read and define body forces
	file.open(init);
	if (file.fail()){
		Ni = 0;
		initial_condition = NULL;
	}
	else
        {
	    file >> Ni; //THE FIRST NUMBER SHOULD BE THE NUMBER OF SOURCES
	    initial_condition = new initial[Ni];
	    for (int i = 0; i<Ni; i++){
		    file >> initial_condition[i].segs[0].point1.x;
			file >> initial_condition[i].segs[0].point1.y;

			file >> initial_condition[i].segs[1].point1.x;
			file >> initial_condition[i].segs[1].point1.y;

			file >> initial_condition[i].segs[2].point1.x;
			file >> initial_condition[i].segs[2].point1.y;

			file >> initial_condition[i].segs[3].point1.x;
			file >> initial_condition[i].segs[3].point1.y;

			file >> initial_condition[i].temp;
			file >> initial_condition[i].temp_eq;

			initial_condition[i].segs[0].point2.x = initial_condition[i].segs[1].point1.x;
			initial_condition[i].segs[0].point2.y = initial_condition[i].segs[1].point1.y;

			initial_condition[i].segs[1].point2.x = initial_condition[i].segs[2].point1.x;
			initial_condition[i].segs[1].point2.y = initial_condition[i].segs[2].point1.y;

			initial_condition[i].segs[2].point2.x = initial_condition[i].segs[3].point1.x;
			initial_condition[i].segs[2].point2.y = initial_condition[i].segs[3].point1.y;

			initial_condition[i].segs[3].point2.x = initial_condition[i].segs[0].point1.x;
			initial_condition[i].segs[3].point2.y = initial_condition[i].segs[0].point1.y;

			// calculate area
			initial_condition[i].area = calc_area(initial_condition[i].segs[0].point1,initial_condition[i].segs[1].point1,initial_condition[i].segs[1].point2)
			+calc_area(initial_condition[i].segs[0].point1,initial_condition[i].segs[2].point1,initial_condition[i].segs[2].point2);

			// center of the shape
			initial_condition[i].center.x = (initial_condition[i].segs[0].point1.x + initial_condition[i].segs[1].point1.x
										  + initial_condition[i].segs[2].point1.x + initial_condition[i].segs[3].point1.x)/4;
			initial_condition[i].center.y = (initial_condition[i].segs[0].point1.y + initial_condition[i].segs[1].point1.y
										  + initial_condition[i].segs[2].point1.y + initial_condition[i].segs[3].point1.y)/4;
	    }
	}
	file.close();
	// total number of sources
	Ntot = Np + Nv + Nb + Ni;
	// array to store the energies
	energies = new double[Ntot];
	source_type = new int[Ntot];
	// define energies associated with prescribed boundaries
	for (int i=0; i<Np; i++)
	{
		if (strcmp(local_type,"TRANSIENT")==0) {
			energies[i] = mat->cumul_CV[mat->Nm-1]*(ptr_to_presc->presc_handle[i].length)*
			(ptr_to_presc->presc_handle[i].temp-ptr_to_presc->presc_handle[i].temp_eq)/4*t_max;
		}
		else
		{
		    energies[i] = mat->cumul_CV[mat->Nm-1]*(ptr_to_presc->presc_handle[i].length)*
		    (ptr_to_presc->presc_handle[i].temp-ptr_to_presc->presc_handle[i].temp_eq)/4;
		}
		source_type[i] = 1;
		cout << "energy of " << i << " : " << energies[i] << endl;
	}
	// define energies associated with volumetric heating
	for (int i=Np; i<Np+Nv; i++)
	{
		if (strcmp(local_type,"TRANSIENT")==0) {
			energies[i] = mat->cumul_C[mat->Nm-1]*(sprd_src_array[i-Np].area)*
		    (sprd_src_array[i-Np].temp-sprd_src_array[i-Np].temp_eq)*t_max;
		}
		else{
		    energies[i] = mat->cumul_C[mat->Nm-1]*(sprd_src_array[i-Np].area)*
		    (sprd_src_array[i-Np].temp-sprd_src_array[i-Np].temp_eq);
		}
		source_type[i] = 2;
		cout << "energy of " << i << " : " << energies[i] << endl;
	}
	// define energies associated with body force
	for (int i=Np+Nv; i<Nv+Np+Nb; i++)
	{
		if (strcmp(local_type,"TRANSIENT")==0) {
			energies[i] = mat->cumul_CV[mat->Nm-1]*(bd_frc_array[i-Np-Nv].area)*
		    sqrt(bd_frc_array[i-Np-Nv].vctr.x*bd_frc_array[i-Np-Nv].vctr.x+bd_frc_array[i-Np-Nv].vctr.y*bd_frc_array[i-Np-Nv].vctr.y)/2*t_max;
		}
		else{
		    energies[i] = mat->cumul_CV[mat->Nm-1]*(bd_frc_array[i-Np-Nv].area)*
		    sqrt(bd_frc_array[i-Np-Nv].vctr.x*bd_frc_array[i-Np-Nv].vctr.x+bd_frc_array[i-Np-Nv].vctr.y*bd_frc_array[i-Np-Nv].vctr.y)/2;
		}
		source_type[i] = 3;
		cout << "energy of " << i << " : " << energies[i] << endl;
	}
	if (strcmp(local_type,"TRANSIENT")==0) {
		for (int i=Np+Nv+Nb; i<Ntot; i++)
		{
			energies[i] = mat->cumul_C[mat->Nm-1]*(initial_condition[i-Np-Nv-Nb].area)*
			(initial_condition[i-Np-Nv-Nb].temp-initial_condition[i-Np-Nv-Nb].temp_eq);
			source_type[i] = 4;
			cout << "energy of " << i << " : " << energies[i] << endl;
		}
	}
    // define energies associated with initial condition

	// calculate and store the absolute cumulative energies
	cumul_energies = new double[Ntot];

	if (Ntot == 0) {
		cout << "There does not seem to be any source => abort !" << endl;
		abort();
	}
	cumul_energies[0] = abs(energies[0]);
	for (int i = 1; i<Ntot; i++){
		cumul_energies[i] = cumul_energies[i-1] + abs(energies[i]);
	}
	total_energy = cumul_energies[Ntot-1];
	N_cumul = new double[Ntot];
	cout << "normalized distribution of sources:" << endl;
	for (int i = 0; i<Ntot; i++){
		N_cumul[i] = cumul_energies[i]/cumul_energies[Ntot-1];
		cout << i << " " <<N_cumul[i] << endl;
	}
}


void sources::display_vol(){
        for (int i=0; i<Nv; i++)
	{
		cout << sprd_src_array[i].segs[0].point1.x << " " << sprd_src_array[i].segs[0].point1.y << " "
		<< sprd_src_array[i].segs[0].point2.x << " " << sprd_src_array[i].segs[0].point2.y << " "
		<< sprd_src_array[i].segs[1].point1.x << " " << sprd_src_array[i].segs[1].point1.y << " "
		<< sprd_src_array[i].segs[1].point2.x << " " << sprd_src_array[i].segs[1].point2.y << " "
		<< sprd_src_array[i].segs[2].point1.x << " " << sprd_src_array[i].segs[2].point1.y << " "
		<< sprd_src_array[i].segs[2].point2.x << " " << sprd_src_array[i].segs[2].point2.y << " "
		<< sprd_src_array[i].segs[3].point1.x << " " << sprd_src_array[i].segs[3].point1.y << " "
		<< sprd_src_array[i].segs[3].point2.x << " " << sprd_src_array[i].segs[3].point2.y << " " << endl;
		cout << sprd_src_array[i].area << " " << sprd_src_array[i].center.x << " " <<sprd_src_array[i].center.y << endl;
		cout << sprd_src_array[i].temp <<  " " << sprd_src_array[i].temp_eq << endl;
	}
}

void sources::display_bod(){
        for (int i=0; i<Nb; i++)
	{
		cout << bd_frc_array[i].segs[0].point1.x << " " << bd_frc_array[i].segs[0].point1.y << " "
		<< bd_frc_array[i].segs[0].point2.x << " " << bd_frc_array[i].segs[0].point2.y << " "
		<< bd_frc_array[i].segs[1].point1.x << " " << bd_frc_array[i].segs[1].point1.y << " "
		<< bd_frc_array[i].segs[1].point2.x << " " << bd_frc_array[i].segs[1].point2.y << " "
		<< bd_frc_array[i].segs[2].point1.x << " " << bd_frc_array[i].segs[2].point1.y << " "
		<< bd_frc_array[i].segs[2].point2.x << " " << bd_frc_array[i].segs[2].point2.y << " "
		<< bd_frc_array[i].segs[3].point1.x << " " << bd_frc_array[i].segs[3].point1.y << " "
		<< bd_frc_array[i].segs[3].point2.x << " " << bd_frc_array[i].segs[3].point2.y << " " << endl;
		cout << bd_frc_array[i].area << " " << bd_frc_array[i].center.x << " " <<  bd_frc_array[i].center.y << endl;
		cout << bd_frc_array[i].vctr.x << " " << bd_frc_array[i].vctr.y << endl;
	}
}

void sources::display_init(){
        for (int i=0; i<Ni; i++)
	{
		cout << initial_condition[i].segs[0].point1.x << " " << initial_condition[i].segs[0].point1.y << " "
		<< initial_condition[i].segs[0].point2.x << " " << initial_condition[i].segs[0].point2.y << " "
		<< initial_condition[i].segs[1].point1.x << " " << initial_condition[i].segs[1].point1.y << " "
		<< initial_condition[i].segs[1].point2.x << " " << initial_condition[i].segs[1].point2.y << " "
		<< initial_condition[i].segs[2].point1.x << " " << initial_condition[i].segs[2].point1.y << " "
		<< initial_condition[i].segs[2].point2.x << " " << initial_condition[i].segs[2].point2.y << " "
		<< initial_condition[i].segs[3].point1.x << " " << initial_condition[i].segs[3].point1.y << " "
		<< initial_condition[i].segs[3].point2.x << " " << initial_condition[i].segs[3].point2.y << " " << endl;
		cout << initial_condition[i].area << " " << initial_condition[i].center.x << " " <<  initial_condition[i].center.y << endl;
	}
}

void sources::emit(particle * part, RandomClass * r) // updates properties of part to make it emitted by the sources
{
	int index_s = choose(r, N_cumul, Ntot);
	double R, phi, nx, ny, Vx, Vy; //those will be useful for calculating the velocity
	// start with the particle sign
	if (energies[index_s]<0)
	{
		part->sig = -1;
	}
	else {
		part->sig = 1;
	}


	int index_m;
	// 3 cases, depending on the index
	if (index_s<Np) { //emission from a precribed temperature wall
		// draw position
        part->pt0 = emit_from_segment(ptr_to_presc->presc_handle[index_s], r);
		//draw mode index
		index_m = choose(r, ptr_to_mat->N_cumul_CV, ptr_to_mat->Nm);
		part->V = ptr_to_mat->VG[index_m]; //norm of velocity
		nx = ptr_to_presc->presc_handle[index_s].normal.x;
		ny = ptr_to_presc->presc_handle[index_s].normal.y;
		R = r->randu(); phi = 2*PI*r->randu(); Vx = sqrt(R)*part->V; Vy = sqrt(1-R)*part->V*cos(phi);
		part->Vp0.x = nx*Vx-ny*Vy;
		part->Vp0.y = ny*Vx+nx*Vy;
		// draw initial time
		part->t = t_max*r->randu();
	}
	else {
		if (index_s<Np+Nv)
		{// emission from volumetric source
			//draw position
            part->pt0 = emit_from_quadrilater(sprd_src_array[index_s-Np], r);
			//draw mode index
			index_m = choose(r, ptr_to_mat->N_cumul_C, ptr_to_mat->Nm);
			part->V = ptr_to_mat->VG[index_m];
			R = 2*r->randu()-1;
			phi = 2*PI*r->randu();
			part->Vp0.x = (part->V)*R; // velocity (2D vector)
			part->Vp0.y = (part->V)*sqrt(1-R*R)*cos(phi);
            // draw initial time
			part->t = t_max*r->randu();


		}
		else {
			if (index_s<Np+Nv+Nb) {
				// emission from body force source
				// if we fall in that case, we need to redefine the sign
				R = r->randu();
				nx = bd_frc_array[index_s-Np-Nv].vctr.x/sqrt(bd_frc_array[index_s-Np-Nv].vctr.x*bd_frc_array[index_s-Np-Nv].vctr.x
														 +bd_frc_array[index_s-Np-Nv].vctr.y*bd_frc_array[index_s-Np-Nv].vctr.y);
				ny = bd_frc_array[index_s-Np-Nv].vctr.y/sqrt(bd_frc_array[index_s-Np-Nv].vctr.x*bd_frc_array[index_s-Np-Nv].vctr.x
														 +bd_frc_array[index_s-Np-Nv].vctr.y*bd_frc_array[index_s-Np-Nv].vctr.y);
				if (R<0.5)
				{
					part->sig = 1;
				}
				else {
					part->sig = -1;
				}

				//draw position
                part->pt0 = emit_from_quadrilater(bd_frc_array[index_s-Np-Nv], r);
				//draw mode index
				index_m = choose(r, ptr_to_mat->N_cumul_CV, ptr_to_mat->Nm);
				part->V = ptr_to_mat->VG[index_m];
				R = r->randu();
				phi = 2*PI*r->randu();
				Vx = -part->sig*(part->V)*sqrt(R); // velocity (2D vector); //"-" sign because the body force is opposite to imposed temperature gradient
				Vy = -part->sig*(part->V)*sqrt(1-R)*cos(phi);
				part->Vp0.x = nx*Vx-ny*Vy;
				part->Vp0.y = ny*Vx+nx*Vy;
				// draw initial time
				part->t = t_max*r->randu();
			}
			else {
				if (index_s<Ntot)
				{
					// emission from initial condition
					//draw position
					part->pt0 = emit_from_quadrilater(initial_condition[index_s-Np-Nv-Nb], r);
					//draw mode index
					index_m = choose(r, ptr_to_mat->N_cumul_C, ptr_to_mat->Nm); // Here we consider that the initial distribution is a Bose-Einstein
					part->V = ptr_to_mat->VG[index_m];
					R = 2*r->randu()-1;
					phi = 2*PI*r->randu();
					part->Vp0.x = (part->V)*R; // velocity (2D vector)
					part->Vp0.y = (part->V)*sqrt(1-R*R)*cos(phi);
					part->t = 0; // time is zero if emitted from initial source
				}
				else {
					cout << "In sources::emit, chosen index not within the expected range" << endl;
					abort();
				}
			}
		}
	}
    //assign the rest of the particle properties
	part->mode_index = index_m;
	part->tau = ptr_to_mat->tau[index_m]; //current relaxation time
	part->pol = ptr_to_mat->pol[index_m]; // current polarization
	part->weight = total_energy/Npart;
	part->alive = true; //born to be alive
	part->counter = 0;
}
/*
 class detectors
 {

 }
 */
double fl(double x, double LX, int N) // "modulo" function (for debugging purpose)
{
    return x-floor(x/LX)*LX;
}

// function used for printing the degree of completion in the main function
void print_percent(int N, int Npart)
{
    if ((100*(long long)N) % Npart ==0)
        cout << (100*(long long)N)/Npart << "% completed" << endl;
}

void matlab_write_geometry(reflective_bdrs * ref, prescribed_bdrs * presc, periodic_bdrs * per, sources * src, detector_array_T * Td, detector_array_H * Hd, const char * filename)
{// this function writes a .m file that may be run in matlab to draw the geometry of the problem of interest
    // open the target file
    ofstream file;
	file.open(filename);
	file << "REFdata = zeros(" << ref->N << ",5);" << endl;
	file << "PREdata = zeros(" << presc->N << ",5);" << endl;
	file << "PERdata = zeros(" << per->N << ",6);" << endl;
	file << "BDdata = zeros(" << src->Nb << ",10);" << endl;
	file << "VOLdata = zeros(" << src->Nv << ",9);" << endl;
	file << "INITdata = zeros(" << src->Ni << ",9);" << endl;
	file << "T_detectors = zeros(" << Td->N << ",8);" << endl;
	file << "H_detectors = zeros(" << Hd->N << ",10);" << endl;

    // read the reflective boundaries
    for (int i=0; i<ref->N; i++)
    {
        file << "REFdata(" << i+1 << ",:)= [ " << ref->ref_handle[i].point1.x << " " << ref->ref_handle[i].point1.y  <<
         " " << ref->ref_handle[i].point2.x << " " << ref->ref_handle[i].point2.y  << " " << ref->ref_handle[i].spec << "];" << endl;
    }
    // read prescribed temperature boundaries
    for (int i=0; i<presc->N; i++)
    {
        file << "PREdata(" << i+1 << ",:)= [ " << presc->presc_handle[i].point1.x << " " << presc->presc_handle[i].point1.y  <<
         " " << presc->presc_handle[i].point2.x << " " << presc->presc_handle[i].point2.y  << " " << presc->presc_handle[i].temp-presc->presc_handle[i].temp_eq << "];" << endl;
    }
    // read periodic boundaries
    for (int i=0; i<per->N; i++)
    {
        file << "PERdata(" << i+1 << ",:)= [ " << per->per_handle[i].point1.x << " " << per->per_handle[i].point1.y  <<
         " " << per->per_handle[i].point2.x << " " << per->per_handle[i].point2.y  << " " << per->per_handle[i].vctr.x << " " << per->per_handle[i].vctr.y << "];" << endl;
    }
    // read source terms
    for (int i = 0; i<src->Nb; i++)
    {
        file << "BDdata(" << i+1 << ",:)= [ " << src->bd_frc_array[i].segs[0].point1.x << " " << src->bd_frc_array[i].segs[0].point1.y  << " " <<
                                               src->bd_frc_array[i].segs[1].point1.x << " " << src->bd_frc_array[i].segs[1].point1.y  << " " <<
                                               src->bd_frc_array[i].segs[2].point1.x << " " << src->bd_frc_array[i].segs[2].point1.y  << " " <<
                                               src->bd_frc_array[i].segs[3].point1.x << " " << src->bd_frc_array[i].segs[3].point1.y  << " " <<
                                   src->bd_frc_array[i].vctr.x << " " << src->bd_frc_array[i].vctr.y << "];" << endl;
    }
    for (int i = 0; i<src->Nv; i++)
    {
        file << "VOLdata(" << i+1 << ",:)= [ " << src->sprd_src_array[i].segs[0].point1.x << " " << src->sprd_src_array[i].segs[0].point1.y  << " " <<
                                               src->sprd_src_array[i].segs[1].point1.x << " " << src->sprd_src_array[i].segs[1].point1.y  << " " <<
                                               src->sprd_src_array[i].segs[2].point1.x << " " << src->sprd_src_array[i].segs[2].point1.y  << " " <<
                                               src->sprd_src_array[i].segs[3].point1.x << " " << src->sprd_src_array[i].segs[3].point1.y  << " " <<
                                   src->sprd_src_array[i].temp-src->sprd_src_array[i].temp_eq << "];" << endl;
    }
    for (int i = 0; i<src->Ni; i++)
    {
        file << "INITdata(" << i+1 << ",:)= [ " << src->initial_condition[i].segs[0].point1.x << " " << src->initial_condition[i].segs[0].point1.y  << " " <<
                                               src->initial_condition[i].segs[1].point1.x << " " << src->initial_condition[i].segs[1].point1.y  << " " <<
                                               src->initial_condition[i].segs[2].point1.x << " " << src->initial_condition[i].segs[2].point1.y  << " " <<
                                               src->initial_condition[i].segs[3].point1.x << " " << src->initial_condition[i].segs[3].point1.y  << " " <<
                                   src->initial_condition[i].temp-src->initial_condition[i].temp_eq << "];" << endl;
    }
    for (int i = 0; i<Td->N; i++)
    {
        file << "T_detectors(" << i+1 << ",:)= [ " << Td->t_handle[i].segs[0].point1.x << " " << Td->t_handle[i].segs[0].point1.y  << " " <<
                                               Td->t_handle[i].segs[1].point1.x << " " << Td->t_handle[i].segs[1].point1.y  << " " <<
                                               Td->t_handle[i].segs[2].point1.x << " " << Td->t_handle[i].segs[2].point1.y  << " " <<
                                               Td->t_handle[i].segs[3].point1.x << " " << Td->t_handle[i].segs[3].point1.y  << " " <<
                                   "];" << endl;
    }
    for (int i = 0; i<Hd->N; i++)
    {
        file << "H_detectors(" << i+1 << ",:)= [ " << Hd->h_handle[i].segs[0].point1.x << " " << Hd->h_handle[i].segs[0].point1.y  << " " <<
                                               Hd->h_handle[i].segs[1].point1.x << " " << Hd->h_handle[i].segs[1].point1.y  << " " <<
                                               Hd->h_handle[i].segs[2].point1.x << " " << Hd->h_handle[i].segs[2].point1.y  << " " <<
                                               Hd->h_handle[i].segs[3].point1.x << " " << Hd->h_handle[i].segs[3].point1.y  << " " <<
                                   Hd->h_handle[i].vctr.x << " " << Hd->h_handle[i].vctr.y  << " " << "];" << endl;
    }
    //code for drawing
    file << "figure; hold on;" << endl;
    file << "for i=1:" << src->Nb << endl;
    file << "hb(i) = fill(BDdata(i,1:2:8),BDdata(i,2:2:8),'g','EdgeColor','none');" << endl;
    file << "end" << endl;
    file << "for i=1:" << src->Nv << endl;
    file << "hv(i) = fill(VOLdata(i,1:2:8),VOLdata(i,2:2:8),'r',LineWidth',2);" << endl;
    file << "end" << endl;
    file << "for i=1:" << src->Ni << endl;
    file << "hi(i) = fill(INITdata(i,1:2:8),INITdata(i,2:2:8),'y','EdgeColor','none');" << endl;
    file << "end" << endl;
    file << "for i=1:" << Td->N << endl;
    file << "ht(i) = fill(T_detectors(i,1:2:8),T_detectors(i,2:2:8),'w','FaceColor','none');" << endl;
    file << "set(ht(i),'LineStyle','--');" << endl;
    file << "end" << endl;
    file << "for i=1:" << Hd->N << endl;
    file << "hh(i) = fill(H_detectors(i,1:2:8),H_detectors(i,2:2:8),'w','FaceColor','none');" << endl;
    file << "set(hh(i),'LineStyle','--');" << endl;
    file << "end" << endl;

    file << "for i=1:" << ref->N << endl;
    file << "plot(REFdata(i,1:2:3),REFdata(i,2:2:4),'k','LineWidth',2);" << endl;
    file << "end" << endl;
    file << "for i=1:" << presc->N << endl;
    file << "plot(PREdata(i,1:2:3),PREdata(i,2:2:4),'b','LineWidth',2);" << endl;
    file << "end" << endl;
    file << "for i=1:" << per->N << endl;
    file << "plot(PERdata(i,1:2:3),PERdata(i,2:2:4),'--','LineWidth',2);" << endl;
    file << "end" << endl;


}



int main()
{
    //INITIALIZE INPUTS
	time_t seconds;
	seconds = time(NULL);
	long seed = (long)seconds;

	RandomClass r;
	r.initialize(seed);

	int NPARTICLES, MAXIMUM_COLL;
	double TEQ;
	ifstream file;
	file.open("input.txt");
	if (file.fail()){
		cout << "NO INPUT FILE" << endl;
	}
	else
        {
	        file >> NPARTICLES;
		file >> MAXIMUM_COLL;
		file >> TEQ;
	}
    // READ GEOMETRY FILES
	reflective_bdrs ref("reflective.txt");
	ref.show();
	prescribed_bdrs presc("prescribed.txt");
	presc.show();
	periodic_bdrs per("periodic.txt");
	per.show();

	// READ MATERIAL PROPERTIES
	materials Si("dataSi.txt", TEQ);
	Si.show_all();

	// READ SOURCE FILES
	sources src(&Si, &presc, "volumetric.txt", "body_force.txt", "initial.txt", NPARTICLES);
	src.display_type();
	cout << src.Np << endl;
	cout << src.Nv << endl;
        cout << src.Nb << endl;
	cout << src.Ni << endl;
	src.display_vol();
	src.display_bod();
	src.display_init();

        particle part(MAXIMUM_COLL);

	// DECLARE AND READ THE HEAT FLUX AND TEMPERATURE DETECTORS
	detector_array_H H_detect("H_detectors.txt", "times.txt");
	detector_array_T T_detect("T_detectors.txt", "times.txt");
	matlab_write_geometry(&ref, &presc, &per, &src, &T_detect, &H_detect, "geometry.m");
	H_detect.show();
	cout << endl;
	T_detect.show();

	// LOOP OVER PARTICLES
	for (int i=0; i<src.Npart; i++){
                print_percent(i, src.Npart);

		src.emit(&part,&r);
  //     cout << part.t << " " << src.t_max << endl;
		while (part.alive){
                        part.initiate_move(&r);

                        part.move(&ref, &presc, &per);

			H_detect.measure(&part);
			T_detect.measure(&part, &Si);

                        part.finish_move(&Si,&r, &ref, &presc, &per);

		}

	}
    // DISPLAY RESULTS IN TERMINAL AND RECORD THEM
	H_detect.show_results();
	H_detect.write("results_H.txt");
	T_detect.show_results();
	T_detect.write("results_T.txt");

	return 0;


}
