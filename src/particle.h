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

#ifndef PARTICLE_H
#define PARTICLE_H

#include "utils.h"
#include "segment.h"
#include "RandomClass.h"
#include "quadrilater.h"
#include "reflective_bdrs.h"
#include "periodic_bdrs.h"
#include "prescribed_bdrs.h"
#include "materials.h"

int choose(RandomClass * r, double norm_cumul[],int N_interv);
extern const double PI;

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
         void move(reflective_bdrs * ref, 
		   prescribed_bdrs * presc, 
		   periodic_bdrs * per); 
	 // uses the information in seg to locate the target point after considering boundaries
        // UPDATES seg (needed for detectors)

         void finish_move(materials * mat, 
			  RandomClass * r,
			  reflective_bdrs * ref, 
			  prescribed_bdrs * presc, 
			  periodic_bdrs * per);
	 void show_all();
         bool P_inside_quad(quadrilater * quad);
 };

#endif
