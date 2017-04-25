/*
Copyright (c) 2017, Jean-Philippe M. Péraud
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

#ifndef UTILS_H
#define UTILS_H

#include "RandomClass.h"

class segment;

struct point{
    double x;
    double y;
};

double det(double ax, double ay, double bx, double by);

double calc_area(point pt1, point pt2, point pt3);

double dist_pts(point pt1, point pt2);

point emit_from_triangle(point pt1, point pt2, point pt3,RandomClass * r);

int choose(RandomClass * r, double norm_cumul[],int N_interv);

double det(segment seg1, segment seg2);

point emit_from_segment(segment seg,RandomClass * r);

int intrsct_test(segment segmentab, segment segment12);
//tests whether two ORIENTED segments intersect
//important: segmentab is a trajectory segment; segment12 is a boundary (any type);
//return +1 if the collision is "direct" (determinant of vectors>0), -1 if indirect, 0 if no collision

point intrsct_pt(segment segmentab, segment segment12);
//calculates the intersection point between two ORIENTED segments intersect
//important: needs to be used along with intrsct_pt to check if the intersection exists

double dist_to_intrsct(segment segmentab,point pt); // calculates the distance between starting point of segment and a point (to be used with previous functions above)

// function used for printing the degree of completion in the main function
void print_percent(int N, int Npart);


#endif
