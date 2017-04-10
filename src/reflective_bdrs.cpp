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

#include "reflective.h"
#include <iostream>
#include <fstream>
#include "reflective_bdrs.h"

reflective_bdrs::reflective_bdrs(const char* filename)
{
	ifstream file;
	file.open(filename);
	if (file.fail())
	{
		cout << "no input file for reflective boundaries" << endl;
		N=0;
		ref_handle = NULL;
	}
	else {
		file >> N;
		ref_handle = new reflective[N];
		for (int i=0; i<N; i++){
			file >> ref_handle[i].point1.x;
			file >> ref_handle[i].point1.y;
			file >> ref_handle[i].point2.x;
			file >> ref_handle[i].point2.y;
			file >> ref_handle[i].spec;
			ref_handle[i].initialize();
		}
	}
}


void reflective_bdrs::show()
{
	for (int i = 0; i<N ; i++){
		cout << ref_handle[i].point1.x << " " << ref_handle[i].point1.y << " "
		    << ref_handle[i].point2.x << " " << ref_handle[i].point2.x << " "
		    << ref_handle[i].spec << " " << ref_handle[i].length << endl;
	}
}

reflective_bdrs::~reflective_bdrs()
{
	delete ref_handle;
}
