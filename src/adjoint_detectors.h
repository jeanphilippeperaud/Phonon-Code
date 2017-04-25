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


#ifndef ADJOINT_DETECTORS_H
#define ADJOINT_DETECTORS_H

class prescribed_bdrs;
class volumetric;
class body_force;
class initial;
class particle;
class sources;

class adjoint_detectors{ // this class will read the relevant source files and define the adjoint detectors to be used.
 public:
  adjoint_detectors(){ptr_to_presc = NULL; sprd_src_array = NULL; bd_frc_array = NULL; initial_condition = NULL; msr_times = NULL;
    steady_H = NULL; steady_T = NULL; H_estimates=NULL; T_estimates = NULL;};
  ~adjoint_detectors();
  adjoint_detectors(prescribed_bdrs * presc, const char * volumetric, const char * body_force, const char * initial, const char * filename_time, sources * src);
  // THIS CONSTRUCTOR MAY ONLY BE CALLED AFTER THE adjoint sources HAVE BEEN INITIALIZED
  void measure(particle * part, sources * src);
  void write(const char * filenameT, const char * filenameH);
  void show_results();
  string type;
  prescribed_bdrs * ptr_to_presc;
  volumetric * sprd_src_array; //pointer to array of volumetric sources
  body_force * bd_frc_array; //pointer to array of body force types of sources
  initial * initial_condition;
  int msr_index;
  double ** H_estimates;
  double ** T_estimates;
  double * steady_H;
  double * steady_T;
  int Np;
  int Nv;
  int Nb;
  int Ni;
  int Ntot;
  int Nt;
  int NT;
  int NH;
  double * msr_times;
};

#endif
