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

#ifndef MATERIALS_H
#define MATERIALS_H

class materials
{
public:
        materials(const char* filename, double Teq); // initialize with file of material properties
        ~materials();
        // including frequencies, densities of states, velocities, relaxation times, Delta frequencies, polarizations.
        // IMPORTANT NOTE: the factor 2 of the TA modes must already be included in the corresponding densities of state.
	int Nm; //total number of modes
	double C; //heat capacity
	double * F; //frequencies
	double * SD; // density of states
	double * VG; //group velocities
	double * tau; // relaxation times
	double * Dom; // delta of frequencies
	int * pol;

	double hbar; //reduced planck constant
	double boltz; //Boltzmann constant
	double * de_dT; // derivative of Bose Einstein

	// useful distributions
	double * C_om;
	double * Ctau;
	double * CV;
	// cumulative distributions (for source terms)
	double * cumul_C;
	double * cumul_Ctau;
	double * cumul_CV;
	// normalized cumulative distributions
	double * N_cumul_C;
	double * N_cumul_Ctau;
	double * N_cumul_CV;

	double kth; // thermal conductivity

	void show_all(); // displays all distributions and data

};

#endif
