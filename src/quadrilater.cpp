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

#include "segment.h"
#include "utils.h"
#include "quadrilater.h"

quadrilater::quadrilater()
{
  segs = new segment[4];
}

quadrilater::quadrilater(segment segment1, segment segment2, segment segment3, segment segment4)
{
  segs = new segment[4];
  segment seg_input[4] = {segment1, segment2, segment3, segment4};
  // assign points to constitutive segments
  for (int i = 0; i<4; i++)
    {
      segs[i].point1.x = seg_input[i].point1.x;
      segs[i].point1.y = seg_input[i].point1.y;
      segs[i].point2.x = seg_input[i].point2.x;
      segs[i].point2.y = seg_input[i].point2.y;
    }

  // calculate area
  area = calc_area(segs[0].point1,segs[1].point1,segs[1].point2)
    +calc_area(segs[0].point1,segs[2].point1,segs[2].point2);
  
  // center of the shape
  center.x = (seg_input[0].point1.x + seg_input[1].point1.x+seg_input[2].point1.x + seg_input[3].point1.x)/4;
  center.y = (seg_input[0].point1.y + seg_input[1].point1.y+seg_input[2].point1.y + seg_input[3].point1.y)/4;
  
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

quadrilater::~quadrilater()
{
  delete[] segs;
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

    if (found1) {
      pt1 = intrsct_pt(seg, quad.segs[i]);
    }
    
    found2 = 0;
    i = -1;
    while (!found2 && i<4){
      i++;
      if (intrsct_test(seg, quad.segs[i])==1) // looking for exit point
	found2 = 1;
    }
    
    if (found1) {
      pt2 = intrsct_pt(seg, quad.segs[i]);
    }


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
