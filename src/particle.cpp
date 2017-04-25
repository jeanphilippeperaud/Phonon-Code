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

#include "particle.h"
#include <math.h>
#include "utils.h"
#define _USE_MATH_DEFINES
 

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
    phi = 2*M_PI*r->randu();
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
	R = r->randu(); phi = 2*M_PI*r->randu(); Vx = sqrt(R)*V; Vy = sqrt(1-R)*V*cos(phi);
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

bool particle::P_inside_quad(quadrilater * quad)
{
  return inside_quad(pt0, quad);
}
