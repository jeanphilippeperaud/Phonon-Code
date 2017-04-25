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

#include <iostream>
#include "particle.h"
#include "volumetric.h"
#include "body_force.h"
#include "initial.h"
#include "materials.h"
#include "prescribed_bdrs.h"
#include "RandomClass.h"
#include "sources.h"
#include "utils.h"
#include <fstream>
#include <string.h>

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
  // The use of this constructor means the use of a "forward" simulation technique
    F_or_B = "FORWARD";
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

sources::sources(materials * mat, const char * T_det, const char * H_det, int NN) // constructor
{
	// this constructor will link the detectors to the sources (adjoint)
    Npart = NN;

    F_or_B = "BACKWARD";
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
	// no "flat" detectors.
	ptr_to_presc = NULL;
	Np = 0;
	// read and define volumetric sources, they come from temperature detectors. However, they are defined only if steady state. If transient, then they have to be counter as initial
	file.open(T_det);
	if (file.fail() || strcmp(type.c_str(),"TRANSIENT")==0){
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

                sprd_src_array[i].temp = 1;
                sprd_src_array[i].temp_eq = 0;// we just set default values for those

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
	    else
	    {
	        sprd_src_array = NULL;
	    }
    }
	file.close();
	//read and define body forces, they come from heat flux detectors
	file.open(H_det);
	if (file.fail()){
		Nb = 0;
		bd_frc_array = NULL;
	}
	else
    {
	    file >> Nb; //THE FIRST NUMBER SHOULD BE THE NUMBER OF SOURCES
	    if (Nb > 0)
        {
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
        else
        {
            bd_frc_array = NULL;
        }
	}
	file.close();
    //read and define initial
	file.open(T_det);
	if (file.fail() || strcmp(type.c_str(),"STEADY")==0){
		Ni = 0;
		initial_condition = NULL;
	}
	else
    {
	    file >> Ni; //THE FIRST NUMBER SHOULD BE THE NUMBER OF SOURCES
	    if (Ni > 0)
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

                initial_condition[i].temp = 1; // this have to somehow be normalized in the adjoint case
                initial_condition[i].temp_eq = 0;

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
	// array to store the energies
	energies = new double[Ntot];
	source_type = new int[Ntot];
	// define energies associated with flat detectors. Not supported in this version but keeping this part of code for future release
	for (int i=0; i<Np; i++)
	{
		if (strcmp(local_type,"TRANSIENT")==0) {
			energies[i] = mat->cumul_CV[mat->Nm-1]*(ptr_to_presc->presc_handle[i].length)*
			(ptr_to_presc->presc_handle[i].temp-ptr_to_presc->presc_handle[i].temp_eq)/4;
		}
		else
		{
		    energies[i] = mat->cumul_CV[mat->Nm-1]*(ptr_to_presc->presc_handle[i].length)*
		    (ptr_to_presc->presc_handle[i].temp-ptr_to_presc->presc_handle[i].temp_eq)/4;
		}
		source_type[i] = 1;
		cout << "energy of " << i << " : " << energies[i] << endl;
	}
	// define energies associated with volumetric heating (steady temperature detector)
	for (int i=Np; i<Np+Nv; i++)
	{
		if (strcmp(local_type,"TRANSIENT")==0) {
			energies[i] = mat->cumul_C[mat->Nm-1]*(sprd_src_array[i-Np].area)*
		    (sprd_src_array[i-Np].temp-sprd_src_array[i-Np].temp_eq)*t_max;
			cout << "EXCEPTION: THIS IS ADJOINT CONSTRUCTOR FOR SOURCES. TRANSIENT MODE DETECTED BUT VOLUMETRIC HEATING NOT NULL" << endl;
			abort();
		}
		else{
		    energies[i] = 1;
		}
		source_type[i] = 2;
		cout << "energy of " << i << " : " << energies[i] << endl;
	}
	// define energies associated with body force
	for (int i=Np+Nv; i<Nv+Np+Nb; i++)
	{
		if (strcmp(local_type,"TRANSIENT")==0) {
			energies[i] = mat->cumul_CV[mat->Nm-1]*
		    sqrt(bd_frc_array[i-Np-Nv].vctr.x*bd_frc_array[i-Np-Nv].vctr.x+bd_frc_array[i-Np-Nv].vctr.y*bd_frc_array[i-Np-Nv].vctr.y)/2;
		}
		else{
		    energies[i] = mat->cumul_CV[mat->Nm-1]*
		    sqrt(bd_frc_array[i-Np-Nv].vctr.x*bd_frc_array[i-Np-Nv].vctr.x+bd_frc_array[i-Np-Nv].vctr.y*bd_frc_array[i-Np-Nv].vctr.y)/2;
		}
		// NOTE: in the adjoint case, transient or steady-state lead to same expression for the energy
		source_type[i] = 3;
		cout << "energy of " << i << " : " << energies[i] << endl;
	}
	if (strcmp(local_type,"TRANSIENT")==0) {
		for (int i=Np+Nv+Nb; i<Ntot; i++)
		{
			energies[i] = 1;
			source_type[i] = 4;
			cout << "energy of " << i << " : " << energies[i] << endl;
		}
	}
    // define energies associated with initial condition

	// calculate and store the absolute cumulative energies
	cumul_energies = new double[Ntot];
	remaining_particles = new int[Ntot];
	if (Ntot == 0) {
		cout << "There does not seem to be any detector => no adjoint source => abort !" << endl;
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
		remaining_particles[i] = (int)floor(NN/Ntot); // each adjoint detector is assigned a number of particles such that the sum of the particles roughly equal the target
		cout << i << " " <<N_cumul[i] << endl;
	}
	cout << "Number of particles per adjoint source: " << floor(NN/Ntot) << endl;
	not_empty = true;
}

void sources::display_vol(){
    cout << "Showing volumetric sources coordinates and intensity: " << endl;
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
	part->pt0 = emit_from_quadrilater(&sprd_src_array[index_s-Np], r);
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
	part->pt0 = emit_from_quadrilater(&bd_frc_array[index_s-Np-Nv], r);
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
	    part->pt0 = emit_from_quadrilater(&initial_condition[index_s-Np-Nv-Nb], r);
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



void sources::emit_adjoint(particle * part, RandomClass * r) // updates properties of part to make it emitted by the sources
{

  double R, phi, nx, ny, Vx, Vy; //those will be useful for calculating the velocity

  // first we need to check that the current emitting adjoint sourc is not "empty"
  if (remaining_particles[src_parser]==0)
    {
      if (src_parser < Ntot-1)
	{
	  src_parser = src_parser+1;
	}
      else {
	not_empty = false;
      }
    }
  if (not_empty){
    int index_m;
    part->sig = 1; // will be different than only for heat flux detectors, half of the cases
    // 3 cases, depending on the index
    if (Np>0)
      {
	cout << "There seems to be a flat detector (not supported yet)" << endl;
	abort();
      }

    if (src_parser < Np+Nv)
      {// emission from volumetric source (corresponding to a temperature detector in steady state)
	//draw position
	part->pt0 = emit_from_quadrilater(&sprd_src_array[src_parser-Np], r);
	// draw mode index
	index_m = choose(r, ptr_to_mat->N_cumul_C, ptr_to_mat->Nm);
	part->V = ptr_to_mat->VG[index_m];
	R = 2*r->randu()-1;
	phi = 2*PI*r->randu();
	part->Vp0.x = (part->V)*R; // velocity (2D vector)
	part->Vp0.y = (part->V)*sqrt(1-R*R)*cos(phi);
	// draw initial time (put 0 because not very important in steady state)
	part->t = 0;

      }
    else {
      if (src_parser<Np+Nv+Nb) {
	// emission from body force source
	// if we fall in that case,  the sign can be either positive or negative
	R = r->randu();
	nx = bd_frc_array[src_parser-Np-Nv].vctr.x/sqrt(bd_frc_array[src_parser-Np-Nv].vctr.x*bd_frc_array[src_parser-Np-Nv].vctr.x
							+bd_frc_array[src_parser-Np-Nv].vctr.y*bd_frc_array[src_parser-Np-Nv].vctr.y);
	ny = bd_frc_array[src_parser-Np-Nv].vctr.y/sqrt(bd_frc_array[src_parser-Np-Nv].vctr.x*bd_frc_array[src_parser-Np-Nv].vctr.x
							+bd_frc_array[src_parser-Np-Nv].vctr.y*bd_frc_array[src_parser-Np-Nv].vctr.y);
	if (R<0.5)
	  {
	    part->sig = 1;
	  }
	else {
	  part->sig = -1;
	}
	//draw position
	part->pt0 = emit_from_quadrilater(&bd_frc_array[src_parser-Np-Nv], r);
	//draw mode index
	index_m = choose(r, ptr_to_mat->N_cumul_CV, ptr_to_mat->Nm);
	part->V = ptr_to_mat->VG[index_m];
	R = r->randu();
	phi = 2*PI*r->randu();
	Vx = -part->sig*(part->V)*sqrt(R); // velocity (2D vector); //no "-" sign because adjoint
	Vy = -part->sig*(part->V)*sqrt(1-R)*cos(phi);
	part->Vp0.x = nx*Vx-ny*Vy;
	part->Vp0.y = ny*Vx+nx*Vy;
	// draw initial time (could either be steady state or not so important to put the transient quantity here)
	part->t = 0;// t_max*r->randu();
      }
      else {
	if (src_parser<Ntot)
	  {
	    // emission from initial condition
	    //draw position
	    part->pt0 = emit_from_quadrilater(&initial_condition[src_parser-Np-Nv-Nb], r);
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

    //assign the rest of the particle properties
    part->mode_index = index_m;
    part->tau = ptr_to_mat->tau[index_m]; //current relaxation time
    part->pol = ptr_to_mat->pol[index_m]; // current polarization
    part->weight = energies[src_parser]/floor(Npart/Ntot);
    part->alive = true; //born to be alive
    part->counter = 0;
  }
  // update the source parser and the remaining  particles in each
  if (remaining_particles[src_parser]>0)
    {
      remaining_particles[src_parser] = remaining_particles[src_parser] - 1;
    }


}
