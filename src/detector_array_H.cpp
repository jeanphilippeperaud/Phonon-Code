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

#include "utils.h"
#include "particle.h"
#include "detector_array_H.h"
#include <fstream>
#include <string.h>
#include "quadrilater.h"


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

//            cout << Vx1 << " " << Vy1 << " " << Vz1 << endl;
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
	    if (inside_quad(pt, &h_handle[i])){
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
      if (overlap_quad(part->seg, &h_handle[i])) {
	cntrbt = overlap_length(part->seg, &h_handle[i]); //this just gives a length
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

