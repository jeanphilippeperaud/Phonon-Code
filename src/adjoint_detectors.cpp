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
#include "sources.h"
#include "volumetric.h"
#include "body_force.h"
#include "particle.h"
#include "prescribed_bdrs.h"
#include "initial.h"
#include <fstream>
#include <string.h>
#include "adjoint_detectors.h"

adjoint_detectors::~adjoint_detectors()  //NOTE: DEBUGGING FLAGS LEFT THERE (AS COMMENTS) IN CASE NEEDED LATER
{
  //        cout <<"delete" << endl;
  /*    for (int i =0 ; i<NT; i++){
	cout << i ;
	for (int j = 0; j<Nt; j++){
	cout << " " << T_estimates[i]+j <<endl;
	}
        cout <<endl;
        }*/
  //cout << "ok 0" << endl;
  if (sprd_src_array!=NULL) {
    //          cout << "ok 00 " << sprd_src_array << endl;
    //   for (int i =0 ; i<Nv; i++){
    //        cout << i << " SPRD " << sprd_src_array+i << endl;
    //    }
    delete[] sprd_src_array;
  }
  //	cout << "ok 1" << endl;
  if (bd_frc_array!=NULL) {
    //	    cout << "ok 01" << endl;
    delete[] bd_frc_array;
  }
  //	cout << "ok 2" << endl;
  if (initial_condition!=NULL) {
    //	    cout << "ok 02" << endl;
    delete[] initial_condition;
  }
  //	cout << "ok 3" << endl;
  if (msr_times!=NULL) {
    //	    cout << "ok 03" << endl;
    //           for (int i =0 ; i<Nt; i++){
    //           cout << i << " MSRT " << msr_times+i << endl;
    //	    }
    delete[] msr_times;
  }
  //	cout << "ok 4" << endl;
  if (steady_H!=NULL) {

    delete[] steady_H;
  }
  if (steady_T!=NULL) {

    delete[] steady_T;
  }
  if (H_estimates!=NULL){
    for (int i=0; i<NH; i++) {
      if (H_estimates[i] != NULL)
	{
	  delete[] H_estimates[i];
	}
    }
    delete[] H_estimates;
  }
  if (T_estimates!=NULL){
    for (int i=0; i<NT; i++) {
      if (T_estimates[i] != NULL)
	{
	  delete[] T_estimates[i];
	}
    }
    delete[] T_estimates;
  }
}

adjoint_detectors::adjoint_detectors(prescribed_bdrs * presc, const char * vol, const char * bod, const char * init, const char * filename_time, sources * src)
{

  //determine the type of simulation (transient of steady)
  int Ntimes;
  point pt1, pt2, pt3, pt4;
  ifstream file, filetime;
  file.open(filename_time);
  if (file.fail()) {
    Nt = 0;
    type = "STEADY";
    msr_times = NULL;
  }
  else{
    type = "TRANSIENT";
    file >> Nt;
    msr_times = new double[Nt];
    msr_index = -1;
    for (int i = 0; i<Nt; i++)
      {
	file >> msr_times[i];
      }
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
      if (Nv>0){
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
      else{
	sprd_src_array = NULL;
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
      if (Nb>0){
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
      else{
	bd_frc_array = NULL;
      }
    }
  file.close();
  //read and define initial conditions
  file.open(init);
  if (file.fail()){
    Ni = 0;
    initial_condition = NULL;
  }
  else
    {
      file >> Ni; //THE FIRST NUMBER SHOULD BE THE NUMBER OF SOURCES
      if (Ni>0)
        {
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
      else
        {
	  initial_condition = NULL;
        }
    }
  file.close();
  // total number of sources
  Ntot = Np + Nv + Nb + Ni;


  //Exception test: This is just for debugging/ bug detection purpose
  if (src->Nv>0 && strcmp(type.c_str(),"TRANSIENT") ==0)
    {
      cout << "detected volumetric adjoint sources although the calculation is set to transient => abort " << endl;
      abort();
    }
  if (src->Ni>0 && strcmp(type.c_str(),"STEADY") ==0)
    {
      cout << "detected initial adjoint sources although the calculation is set to steady => abort " << endl;
      abort();
    }
  // Allocate and initialize the estimates from information provided by the adjoint source
  // normally the size of the temperature estimate array is given by the number of volumetric + initial sources
  NT = src->Nv + src->Ni;
  // the size of the heatflux estimates array is given by the number of body force sources
  NH = src->Nb;
  if (NT>0 && Nt>0) // transient case with T detectors
    {
      T_estimates = new double*[NT];
      for (int i = 0; i<NT; i++){
	T_estimates[i] = new double[Nt];
	for (int j=0; j<Nt; j++) {
	  T_estimates[i][j] = 0;
	}
      }
    }
  if (NH>0 && Nt>0) // transient case with H detectors
    {
      H_estimates = new double*[NH];
      for (int i = 0; i<NH; i++){
	H_estimates[i] = new double[Nt];
	for (int j=0; j<Nt; j++) {
	  H_estimates[i][j] = 0;
	}
      }
    }
  if (NH>0 && Nt == 0)
    {
      steady_H = new double[NH];
      for (int i=0; i<NH; i++)
	{
	  steady_H[i]=0;
	}
    }
  if (NT>0 && Nt == 0)
    {
      steady_T = new double[NT];
      for (int i=0; i<NT; i++)
    	{
	  steady_T[i]=0;
	}
    }
}

void adjoint_detectors::measure(particle * part, sources * src) // measure with the adjoint detectors
{
  double cntrbt;
  int type_dect;
  int index_dect;
  //  cout << "flag 1" << endl;
  point pt, nrmlzd_seg; // point structure to store the coordinates of the normalized vector colinear with segment;
  segment seg1;
  // First, we define the type of detector we are currently looking at and the associated index
  if (src->src_parser < src->Nv+src->Np) { // this one means the adjoint source corresponds to volumetric i.e. temperature detector in steady case
    type_dect = 1;
    index_dect = src->src_parser-src->Np;
  }
  if (src->Nv+src->Np <= src->src_parser &&src->src_parser <src->Nv+src->Nb+src->Np) { // this one means heat flux detector
    type_dect = 2;
    index_dect = src->src_parser-(src->Nv+src->Np);
  }
  if (src->src_parser>=src->Ntot-src->Ni) { // this one means temperature detector, transient case
    type_dect = 3;
    index_dect = src->src_parser-(src->Nv+src->Np+src->Nb);
  }
  //    cout << "flag 2" << endl;
  // Second, we check if there is a contribution from contact with emitting wall
  if (part->collision_type == 2)
    {
      // if transient: we need to add the contribution to ALL times below the collision time
      if (strcmp(type.c_str(),"TRANSIENT")==0) {
	for (int j = 0; j<Nt; j++) {
	  if (part->t + part->Dt <msr_times[j]) {
	    if (type_dect == 1) // all cases should be the same but splitting for clarity
	      {
		T_estimates[index_dect][j] = T_estimates[index_dect][j] + part->sig*part->weight*
		  (ptr_to_presc->presc_handle[part->collision_index].temp-ptr_to_presc->presc_handle[part->collision_index].temp_eq);
		cout << "DEBUGGING: this case should not happen 1" << endl;
	      }
	    if (type_dect == 2) {
	      H_estimates[index_dect][j] = H_estimates[index_dect][j] + part->sig*part->weight*
		(ptr_to_presc->presc_handle[part->collision_index].temp-ptr_to_presc->presc_handle[part->collision_index].temp_eq);
	    }
	    if (type_dect == 3) {
	      T_estimates[index_dect][j] = T_estimates[index_dect][j] + part->sig*part->weight*
		(ptr_to_presc->presc_handle[part->collision_index].temp-ptr_to_presc->presc_handle[part->collision_index].temp_eq);

	    }
	  }
	}
      }
      else {
	// if steady: just need to add the contribution
	if (type_dect == 1) //steady case
	  {
	    steady_T[index_dect] = steady_T[index_dect] + part->sig*part->weight*
	      (ptr_to_presc->presc_handle[part->collision_index].temp-ptr_to_presc->presc_handle[part->collision_index].temp_eq);
	  }
	if (type_dect == 2) {
	  steady_H[index_dect] = steady_H[index_dect] + part->sig*part->weight*
	    (ptr_to_presc->presc_handle[part->collision_index].temp-ptr_to_presc->presc_handle[part->collision_index].temp_eq);
	}
	if (type_dect == 3) {
	  steady_T[index_dect] = steady_T[index_dect] + part->sig*part->weight*
	    (ptr_to_presc->presc_handle[part->collision_index].temp-ptr_to_presc->presc_handle[part->collision_index].temp_eq);
	  cout << "DEBUGGING: this case should not happen 2" << endl;
	}
      }
    }

  //cout << "flag 3" << endl;
  for (int i  = 0; i<Nv; i++) //run through all volumetric sources (used as adjoint detectors)
    {
      if (strcmp(type.c_str(),"TRANSIENT")==0) {
	seg1.point1 = part->pt0;
	for (int j = 0; j<Nt; j++) {
	  if (part->t < msr_times[j] && part->t+part->Dt > msr_times[j]) { // for those times, it will only count partially
	    // define point at given times
	    pt.x = part->pt0.x + (msr_times[j] - part->t)*part->Vp0.x;
	    pt.y = part->pt0.y + (msr_times[j] - part->t)*part->Vp0.y;
	    seg1.point2 = pt;
	    //calculate the contribution based on the overlap between the segment and the adjoint detector
	    // cout << "flag 32" << endl;
	    if (overlap_quad(seg1, &sprd_src_array[i])) // first, test if there is any overlap
	      {
		//   cout << "flag 33" << endl;
		cntrbt = overlap_length(seg1, &sprd_src_array[i]); // this gives a (2D) distance
		if (type_dect == 1) // all cases should be the same but splitting for clarity
		  {
		    T_estimates[index_dect][j] = T_estimates[index_dect][j] + part->sig*part->weight*cntrbt/sqrt(part->Vp0.x*part->Vp0.x+part->Vp0.y*part->Vp0.y)*
		      (sprd_src_array[i].temp-sprd_src_array[i].temp_eq);
		    // unit is temperature
		    cout << "DEBUGGING: this case should not happen 3" << endl;
		  }
		if (type_dect == 2) {
		  H_estimates[index_dect][j] = H_estimates[index_dect][j] + part->sig*part->weight*cntrbt/sqrt(part->Vp0.x*part->Vp0.x+part->Vp0.y*part->Vp0.y)*
		    (sprd_src_array[i].temp-sprd_src_array[i].temp_eq);
		  // unit is W/m^2
		}
		if (type_dect == 3) {
		  T_estimates[index_dect][j] = T_estimates[index_dect][j] + part->sig*part->weight*cntrbt/sqrt(part->Vp0.x*part->Vp0.x+part->Vp0.y*part->Vp0.y)*
		    (sprd_src_array[i].temp-sprd_src_array[i].temp_eq);
		}
	      }
	  }
	  if (part->t+part->Dt < msr_times[j]) { // for those times it will count entirely: the "ending" point is just part->pt1
	    seg1.point2 = part->pt1;
	    //calculate the contribution based on the overlap between the segment and the adjoint detector
	    if (overlap_quad(seg1, &sprd_src_array[i])) // first, test if there is any overlap
	      {
		cntrbt = overlap_length(seg1, &sprd_src_array[i]); // this gives a (2D) distance
		if (type_dect == 1) // all cases should be the same but splitting for clarity
		  {
		    T_estimates[index_dect][j] = T_estimates[index_dect][j] + part->sig*part->weight*cntrbt/sqrt(part->Vp0.x*part->Vp0.x+part->Vp0.y*part->Vp0.y)*
		      (sprd_src_array[i].temp-sprd_src_array[i].temp_eq);
		    // unit is temperature
		    cout << "DEBUGGING: this case should not happen 3" << endl;
		  }
		if (type_dect == 2) {
		  H_estimates[index_dect][j] = H_estimates[index_dect][j] + part->sig*part->weight*cntrbt/sqrt(part->Vp0.x*part->Vp0.x+part->Vp0.y*part->Vp0.y)*
		    (sprd_src_array[i].temp-sprd_src_array[i].temp_eq);
		  // unit is W/m^2
		}
		if (type_dect == 3) {
		  T_estimates[index_dect][j] = T_estimates[index_dect][j] + part->sig*part->weight*cntrbt/sqrt(part->Vp0.x*part->Vp0.x+part->Vp0.y*part->Vp0.y)*
		    (sprd_src_array[i].temp-sprd_src_array[i].temp_eq);
		}
	      }
	  }
	}
      }
      else {
	// if steady: just need to add the contribution
	seg1.point1 = part->pt0;
	seg1.point2 = part->pt1;

	//calculate the contribution based on the overlap between the segment and the adjoint detector
	if (overlap_quad(seg1, &sprd_src_array[i])) // first, test if there is any overlap
	  {
	    cntrbt = overlap_length(seg1, &sprd_src_array[i]); // this gives a (2D) distance
	    if (type_dect == 1) //steady case
	      {
		steady_T[index_dect] = steady_T[index_dect] + part->sig*part->weight*cntrbt/sqrt(part->Vp0.x*part->Vp0.x+part->Vp0.y*part->Vp0.y)*
		  (sprd_src_array[i].temp-sprd_src_array[i].temp_eq);
	      }
	    if (type_dect == 2) {
	      steady_H[index_dect] = steady_H[index_dect] + part->sig*part->weight*cntrbt/sqrt(part->Vp0.x*part->Vp0.x+part->Vp0.y*part->Vp0.y)*
		(sprd_src_array[i].temp-sprd_src_array[i].temp_eq);
	    }
	    if (type_dect == 3) {
	      steady_T[index_dect] = steady_T[index_dect] + part->sig*part->weight*cntrbt/sqrt(part->Vp0.x*part->Vp0.x+part->Vp0.y*part->Vp0.y)*
		(sprd_src_array[i].temp-sprd_src_array[i].temp_eq);
	      cout << "DEBUGGING: this case should not happen 4" << endl;
	    }
	  }
      }
    }

  //cout << "flag 4" << endl;

  for (int i  = 0; i<Nb; i++) //run through all body force sources (used as adjoint detectors)
    {
      if (strcmp(type.c_str(),"TRANSIENT")==0) {

	for (int j = 0; j<Nt; j++) {
	  if (part->t <= msr_times[j] && part->t+part->Dt > msr_times[j]) { // for those times, it will only count partially
	    // define point at given times
	    pt.x = part->pt0.x + (msr_times[j] - part->t)*part->Vp0.x;
	    pt.y = part->pt0.y + (msr_times[j] - part->t)*part->Vp0.y;
	    seg1.point1 = part->pt0;
	    seg1.point2 = pt;
	    //calculate the contribution based on the overlap between the segment and the adjoint detector
	    if (overlap_quad(seg1, &bd_frc_array[i])) // first, test if there is any overlap
	      {
		cntrbt = overlap_length(seg1, &bd_frc_array[i]); // this gives a (2D) distance
		nrmlzd_seg.x = (part->seg.point2.x - part->seg.point1.x)/
		  sqrt((part->seg.point2.x - part->seg.point1.x)*(part->seg.point2.x - part->seg.point1.x)+(part->seg.point2.y - part->seg.point1.y)*(part->seg.point2.y - part->seg.point1.y));
		nrmlzd_seg.y = (part->seg.point2.y - part->seg.point1.y)/
		  sqrt((part->seg.point2.x - part->seg.point1.x)*(part->seg.point2.x - part->seg.point1.x)+(part->seg.point2.y - part->seg.point1.y)*(part->seg.point2.y - part->seg.point1.y));

		cntrbt = cntrbt*(nrmlzd_seg.x*bd_frc_array[i].vctr.x + nrmlzd_seg.y*bd_frc_array[i].vctr.y);

		if (type_dect == 1) // all cases should be the same but splitting for clarity
		  {
		    T_estimates[index_dect][j] = T_estimates[index_dect][j] + part->sig*part->weight*cntrbt;
		    // unit is temperature
		    cout << "DEBUGGING: this case should not happen 3" << endl;
		  }
		if (type_dect == 2) {
		  H_estimates[index_dect][j] = H_estimates[index_dect][j] + part->sig*part->weight*cntrbt;
		  // unit is W/m^2
		}
		if (type_dect == 3) {
		  T_estimates[index_dect][j] = T_estimates[index_dect][j] + part->sig*part->weight*cntrbt;
		}
	      }
	  }
	  if (part->t+part->Dt <= msr_times[j]) { // for those times it will count entirely: the "ending" point is just part->pt1
	    seg1.point1 = part->pt0;
	    seg1.point2 = part->pt1;

	    //calculate the contribution based on the overlap between the segment and the adjoint detector
	    if (overlap_quad(seg1, &bd_frc_array[i])) // first, test if there is any overlap
	      {
		cntrbt = overlap_length(seg1, &bd_frc_array[i]); // this gives a (2D) distance
		nrmlzd_seg.x = (part->seg.point2.x - part->seg.point1.x)/
		  sqrt((part->seg.point2.x - part->seg.point1.x)*(part->seg.point2.x - part->seg.point1.x)+(part->seg.point2.y - part->seg.point1.y)*(part->seg.point2.y - part->seg.point1.y));
		nrmlzd_seg.y = (part->seg.point2.y - part->seg.point1.y)/
		  sqrt((part->seg.point2.x - part->seg.point1.x)*(part->seg.point2.x - part->seg.point1.x)+(part->seg.point2.y - part->seg.point1.y)*(part->seg.point2.y - part->seg.point1.y));
		cntrbt = cntrbt*(nrmlzd_seg.x*bd_frc_array[i].vctr.x + nrmlzd_seg.y*bd_frc_array[i].vctr.y);

		if (type_dect == 1) // all cases should be the same but splitting for clarity
		  {
		    T_estimates[index_dect][j] = T_estimates[index_dect][j] + part->sig*part->weight*cntrbt;
		    // unit is temperature
		    cout << "DEBUGGING: this case should not happen 5" << endl;
		  }
		if (type_dect == 2) {
		  H_estimates[index_dect][j] = H_estimates[index_dect][j] + part->sig*part->weight*cntrbt;
		  // unit is W/m^2
		}
		if (type_dect == 3) {
		  T_estimates[index_dect][j] = T_estimates[index_dect][j] + part->sig*part->weight*cntrbt;
		}
	      }
	  }
	}
      }
      else {
	// if steady: just need to add the contribution
	seg1.point1 = part->pt0;
	seg1.point2 = part->pt1;

	//calculate the contribution based on the overlap between the segment and the adjoint detector
	if (overlap_quad(seg1, &bd_frc_array[i])) // first, test if there is any overlap
	  {
	    cntrbt = overlap_length(seg1, &bd_frc_array[i]);
	    nrmlzd_seg.x = (part->seg.point2.x - part->seg.point1.x)/
	      sqrt((part->seg.point2.x - part->seg.point1.x)*(part->seg.point2.x - part->seg.point1.x)+(part->seg.point2.y - part->seg.point1.y)*(part->seg.point2.y - part->seg.point1.y));
	    nrmlzd_seg.y = (part->seg.point2.y - part->seg.point1.y)/
	      sqrt((part->seg.point2.x - part->seg.point1.x)*(part->seg.point2.x - part->seg.point1.x)+(part->seg.point2.y - part->seg.point1.y)*(part->seg.point2.y - part->seg.point1.y));
	    cntrbt = cntrbt*(nrmlzd_seg.x*bd_frc_array[i].vctr.x + nrmlzd_seg.y*bd_frc_array[i].vctr.y);
	    if (type_dect == 1) //steady case
	      {
		steady_T[index_dect] = steady_T[index_dect] + part->sig*part->weight*cntrbt;
	      }
	    if (type_dect == 2) {
	      steady_H[index_dect] = steady_H[index_dect] + part->sig*part->weight*cntrbt;
	    }
	    if (type_dect == 3) {
	      steady_T[index_dect] = steady_T[index_dect] + part->sig*part->weight*cntrbt;
	      cout << "DEBUGGING: this case should not happen 6" << endl;
	    }
	  }
      }
    }


  //cout << "flag 5" << endl;

  for (int i = 0; i<Ni; i++) //run through all INITIAL CONDITIONS
    {
      //NOTE: if any of those exist, this means this is a transient calculation
      if (strcmp(type.c_str(),"STEADY")==0){
	cout << "Warning: detected data defining an initial condition although the simulation mode is STEADY => ignoring the initial condition" << endl;
	cout << "Relaunch without initial condition to get rid of this annoying message" << endl;
      }
      else{
	for (int j = 0; j<Nt; j++){
	  if (msr_times[j] >=part->t && msr_times[j] < part->t+part->Dt){ // simply check the measurement times between the particle last and next collision
	    pt.x = part->pt0.x + part->Vp0.x*(msr_times[j]-part->t);
	    pt.y = part->pt0.y + part->Vp0.y*(msr_times[j]-part->t);
	    // check whether the particle would be in any of the detectors at these moments

	    if (inside_quad(pt, &initial_condition[i])){
	      if (type_dect == 1) // all cases should be the same but splitting for clarity
		{
		  T_estimates[index_dect][j] = T_estimates[index_dect][j] + part->sig*part->weight*(initial_condition[i].temp-initial_condition[i].temp_eq);
		  // unit is temperature
		  cout << "DEBUGGING: this case should not happen 7" << endl;
		}
	      if (type_dect == 2) {
		H_estimates[index_dect][j] = H_estimates[index_dect][j] + part->sig*part->weight*(initial_condition[i].temp-initial_condition[i].temp_eq);
		// unit is W/m^2
	      }
	      if (type_dect == 3) {
		T_estimates[index_dect][j] = T_estimates[index_dect][j] + part->sig*part->weight*(initial_condition[i].temp-initial_condition[i].temp_eq);
	      }

	    }
	  }

	}
      }
    }

  if (strcmp(type.c_str(), "TRANSIENT")==0 && part->t + part->Dt > msr_times[Nt-1]) // the measure funtion handles this termination case because the detector knows the measuring times.
    {
      msr_index = -1; // only for transient cases: reset the msr_index
      part->alive = 0; // "kill" the particle
    }

}

void adjoint_detectors::write(const char * filenameT, const char * filenameH) //the first argument is the name where the temperature estimates are written
// the second argument is the name where heat flux estimates are written.
{
  ofstream file;
  file.open(filenameT);
  for (int i = 0; i<NT ; i++){
    if (strcmp(type.c_str(),"TRANSIENT")==0){
      for (int j = 0; j<Nt; j++)
	{
	  file << T_estimates[i][j];
	  file << " " ;
	}
      file << endl;
    }
    else{
      file << steady_T[i] << endl;
    }
  }
  file.close();

  file.open(filenameH);
  for (int i = 0; i<NH; i++)
    {
      if (strcmp(type.c_str(),"TRANSIENT")==0){
	for (int j = 0; j<Nt; j++)
	  {
	    file << H_estimates[i][j];
	    file << " " ;
	  }
	file << endl;
      }
      else{
	file << steady_H[i] << endl;
      }
    }
}

void adjoint_detectors::show_results() //displays the results
{
  if (NT>0)
    cout << "Showing temperature estimates: " << endl;
  for (int i = 0; i<NT ; i++){
    if (strcmp(type.c_str(),"TRANSIENT")==0){
      for (int j = 0; j<Nt; j++)
	{
	  cout << T_estimates[i][j];
	  cout << " " ;
	}
      cout << endl;
    }
    else{
      cout << steady_T[i] << endl;
    }
  }
  cout << endl;
  if (NH>0)
    cout << "Showing heat flux estimates: " << endl;
  for (int i = 0; i<NH; i++)
    {
      if (strcmp(type.c_str(),"TRANSIENT")==0){
	for (int j = 0; j<Nt; j++)
	  {
	    cout << H_estimates[i][j];
	    cout << " " ;
	  }
	cout << endl;
      }
      else{
	cout << steady_H[i] << endl;
      }
    }
  cout << endl;
}
