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



#ifndef SOURCES_H
#define SOURCES_H

extern const double PI;

class materials;
class prescribed_bdrs;
class volumetric;
class body_force;
class initial;
class particle;
class RandomClass;

class sources
{
 public:
  sources(materials * mat, prescribed_bdrs * presc, const char * vol, const char * bod, const char * init, int NN); // vol and bod are respectively the name files of volumetric sources and body force sources
  sources(materials * mat, const char * T_det, const char * H_det, int NN); // constructor used for adjoint simulations
  ~sources();
  string type; //TRANSIENT or STEADY
  string F_or_B;
  double total_energy;
  double t_max;  // maximum time (for TRANSIENT only)
  double * energies; // individual weights of each source
  int * remaining_particles; // for the ADJOINT version, a counter is implemented for each adjoint source. Exactly Npart particles are emitted from EACH source
  double * cumul_energies; // same, but cumulative
  double * N_cumul; // same, but normalized
  int * source_type; // for determining from which source to emit 1 = prescribed temperature wall; 2 = volumetric heating; 3 = body force induced by temperature gradient
  prescribed_bdrs * ptr_to_presc;
  volumetric * sprd_src_array; //pointer to array of volumetric sources
  body_force * bd_frc_array; //pointer to array of body force types of sources
  initial * initial_condition;
  int src_parser; //indicates which adjoint source we are currently treating
  bool not_empty; // switches to false when src_parser reaches the number of adjoint sources and the corresponding remaining particles become zero
  int Np; //number of prescribed sources
  int Nv; // number of volumetric sources
  int Nb; // number of body_force sources
  int Ni; // number of sources associated with initial condition
  int Ntot; //total number of separate sources
  int Npart; // total number of particles to be simulated
  materials * ptr_to_mat;
  void emit(particle * part, RandomClass * r); // returns a particle sampled from the sources
  void emit_adjoint(particle * part, RandomClass * r); // return a particle sampled from the adjoint sources (detectors). This should only be used if the
// proper constructor (overloaded for adjoint simulations) has been called
  void display_type(){cout << type << endl;} // to check if this is transient or steady simulation
  void display_vol();
  void display_bod();
  void display_init();
};

#endif
