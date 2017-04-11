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
#include <fstream>
#include <math.h>
#include "materials.h"

materials::materials(const char* filename, double Teq)
{ // filename must obey very strict data organization, as follows:
// 1st entry = total number of "modes" (or "frequency bins")
// then each line = frequency, density of states, velocity, Delta frequency, relaxation time, polarization
	// Teq is the equilibrium temperature, necessary for calculating the temperature dependent properties
	hbar = 1.054517e-34;
	boltz = 1.38065e-23;
	std::ifstream file;
	file.open(filename);
	if (file.fail())
	{
		std::cout << "no input file for material data" << std::endl;
		Nm=0;
		F = NULL;
		SD = NULL;
		VG = NULL;
		tau = NULL;
		Dom = NULL;
		pol = NULL;
		de_dT = NULL;
		C_om = NULL;
		Ctau = NULL;
		CV = NULL;
		cumul_C = NULL;
		cumul_Ctau = NULL;
		cumul_CV = NULL;
		N_cumul_C = NULL;
		N_cumul_Ctau = NULL;
		N_cumul_CV = NULL;
	}
	else {
		file >> Nm;
		F = new double[Nm];
		SD = new double[Nm];
		VG = new double[Nm];
		tau = new double[Nm];
		Dom = new double[Nm];
		pol = new int[Nm];
		de_dT = new double[Nm];
		C_om = new double[Nm];
		Ctau = new double[Nm];
		CV = new double[Nm];
		cumul_C = new double[Nm];
		cumul_Ctau = new double[Nm];
		cumul_CV = new double[Nm];
		N_cumul_C = new double[Nm];
		N_cumul_Ctau = new double[Nm];
		N_cumul_CV = new double[Nm];
		// read material parameters
		for (int i=0; i<Nm; i++){
			file >> F[i];
			file >> SD[i];
			file >> VG[i];
			file >> Dom[i];
			file >> tau[i];
			file >> pol[i];
			de_dT[i] = (hbar*F[i]/Teq)*(hbar*F[i]/Teq)/boltz*exp(hbar*F[i]/boltz/Teq)/(exp(hbar*F[i]/boltz/Teq)-1)/(exp(hbar*F[i]/boltz/Teq)-1);
			C_om[i] = SD[i]*de_dT[i]*Dom[i];
			Ctau[i] = SD[i]*de_dT[i]*Dom[i]/tau[i];
			CV[i] = SD[i]*de_dT[i]*Dom[i]*VG[i];
		}

		// calculate cumulative distributions
		cumul_C[0] = C_om[0];
		cumul_Ctau[0] = Ctau[0];
		cumul_CV[0] = CV[0];
		kth = CV[0]*VG[0]*tau[0]/3;
		for (int i=1; i<Nm; i++){
			cumul_C[i] = cumul_C[i-1] + C_om[i];
			cumul_Ctau[i] = cumul_Ctau[i-1] + Ctau[i];
			cumul_CV[i] = cumul_CV[i-1] + CV[i];
			kth = kth + CV[i]*VG[i]*tau[i]/3;
		}

		// normalize cumulative distributions
		for (int i = 0; i<Nm; i++)
		{
			N_cumul_C[i] = cumul_C[i]/cumul_C[Nm-1];
			N_cumul_Ctau[i] = cumul_Ctau[i]/cumul_Ctau[Nm-1];
			N_cumul_CV[i] = cumul_CV[i]/cumul_CV[Nm-1];
		}
		C = cumul_C[Nm-1];
	}
}

materials::~materials()
{
	delete[] F;
	delete[] SD;
	delete[] VG;
	delete[] tau;
	delete[] Dom;
	delete[] pol;
	delete[] de_dT;
	delete[] C_om;
	delete[] Ctau;
	delete[] CV;
	delete[] cumul_C;
	delete[] cumul_Ctau;
	delete[] cumul_CV;
	delete[] N_cumul_C;
	delete[] N_cumul_Ctau;
	delete[] N_cumul_CV;
}

void materials::show_all()
{
        std::cout << "number of frequency cells: " << Nm << std::endl;
	std::cout << "Boltzmann constant: " << boltz << std::endl;
	std::cout << "reduced Planck constant: " << hbar << std::endl;
	std::cout << "RAW DATA: " << std::endl;
	for (int i = 0; i<Nm; i++){
		std::cout << F[i] << " " << SD[i] << " " << VG[i] << " " << Dom[i] << " " << tau[i] << " " << pol[i] << std::endl;
	}
	std::cout << "DISTRIBUTIONS: " << std::endl;
	for (int i = 0; i<Nm; i++){
		std::cout << i << " " <<C_om[i] << " " << Ctau[i] << " " << CV[i] << std::endl;
	}
	std::cout << "CUMULATIVE DISTRIBUTIONS: " << std::endl;
	for (int i = 0; i<Nm; i++){
		std::cout << i << " " <<cumul_C[i] << " " << cumul_Ctau[i] << " " << cumul_CV[i] << std::endl;
	}
	std::cout << "NORMALIZED CUMULATIVE DISTRIBUTIONS:" << std::endl;
	for (int i = 0; i<Nm; i++){
		std::cout << i << " " << N_cumul_C[i] << " " << N_cumul_Ctau[i] << " " << N_cumul_CV[i] << std::endl;
	}
	std::cout << "thermal conductivity : " << kth << std::endl;

}


