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

#include "detector_T.h"
#include "utils.h"
#include "segment.h"

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

