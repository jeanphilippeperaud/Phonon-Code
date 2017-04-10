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

#ifndef QUADRILATER_H
#define QUADRILATER_H
#include "segment.h"
#include "utils.h"

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


int intrsct_test(segment segmentab, segment segment12);
point intrsct_pt(segment segmentab, segment segment12);

bool inside_quad(point pt, quadrilater quad); // checks if point pt is inside quad. IMPORTANT: only works if quad is convex

bool overlap_quad(segment seg, quadrilater quad);
// tests whether a segment and a quadrilater overlap

double overlap_length(segment seg, quadrilater quad); // calculates the length of segment overlap
// IMPORTANT: if this function is used although there is no overlap, it can return a non zero value. The function overlap_quad
// MUST be used before using this function.

point emit_from_quadrilater(quadrilater quad, RandomClass * r);
// subdivide into 2 triangles and calculate their respective absolute areas


#endif
