#include "utils.h"
#include "segment.h"
#include "RandomClass.h"
#include "reflective_bdrs.h"
#include "prescribed_bdrs.h"
#include "periodic_bdrs.h"
#include "sources.h"
#include "materials.h"
#include "particle.h"
#include "adjoint_detectors.h"
#include "detector_array_T.h"
#include "detector_array_H.h"
#include "plot_utils.h"
#include <fstream>
#include <string.h>
#include <stdio.h>

const double PI=3.141592653589;

int main()
{    //INITIALIZE INPUTS
  time_t seconds;
  seconds = time(NULL);
  long seed = (long)seconds;

  RandomClass r;
  r.initialize(seed);

  string SIMU_TYPE;
  int NPARTICLES, MAXIMUM_COLL;
  double TEQ;
  ifstream file;
  file.open("input.txt");
  if (file.fail()){
    cout << "NO INPUT FILE" << endl;
  }
  else
    {
      file >> SIMU_TYPE;  //case sensitive
      file >> NPARTICLES;
      file >> MAXIMUM_COLL;
      file >> TEQ;
    }

  // READ GEOMETRY FILES
  reflective_bdrs ref("reflective.txt");
  ref.show();
  prescribed_bdrs presc("prescribed.txt");
  presc.show();
  periodic_bdrs per("periodic.txt");
  per.show();

  // READ MATERIAL PROPERTIES
  materials Si("dataSi.txt", TEQ);
  Si.show_all();

  // DECLARE AND READ THE HEAT FLUX AND TEMPERATURE DETECTORS
  detector_array_H H_detect("H_detectors.txt", "times.txt");
  detector_array_T T_detect("T_detectors.txt", "times.txt");

  // PRINT GEOMETRY IN FILE
  particle part(MAXIMUM_COLL);

  // READ SOURCE FILES
  sources src(&Si, &presc, "volumetric.txt", "body_force.txt", "initial.txt", NPARTICLES);
  src.display_type();
  cout << src.Np << endl;
  cout << src.Nv << endl;
  cout << src.Nb << endl;
  cout << src.Ni << endl;
  src.display_vol();
  src.display_bod();
  src.display_init();

  matlab_write_geometry(&ref, &presc, &per, &src, &T_detect, &H_detect, "geometry.m");

  if (strcmp(SIMU_TYPE.c_str(),"BOTH")==0 || strcmp(SIMU_TYPE.c_str(),"BACKWARD")==0 || strcmp(SIMU_TYPE.c_str(),"ADJOINT")==0)
    {

      // READ SOURCE FILES
      sources src_adj(&Si, "T_detectors.txt", "H_detectors.txt", NPARTICLES);

      adjoint_detectors adj_det(&presc, "volumetric.txt", "body_force.txt", "initial.txt", "times.txt", &src_adj);
      src.display_type();

      H_detect.show();
      cout << endl;
      T_detect.show();

      cout << endl << "STARTING ADJOINT SIMULATION" << endl;

      // LOOP OVER PARTICLES
      cout << src_adj.Nv << endl;
      while (src_adj.not_empty) {
	if (src_adj.remaining_particles[src_adj.src_parser]==1){
	  cout << src_adj.src_parser + 1 << " out of " << src_adj.Ntot << " adjoint sources " << endl;
	}


	src_adj.emit_adjoint(&part,&r);

	while (part.alive){
	  part.initiate_move(&r);

	  part.move(&ref, &presc, &per);

	  //H_detect.measure(&part);
	  //T_detect.measure(&part, &Si);
	  adj_det.measure(&part, &src_adj);
	  part.finish_move(&Si,&r, &ref, &presc, &per);

	}

      }
      // DISPLAY RESULTS IN TERMINAL AND RECORD THEM
      adj_det.show_results();
      adj_det.write("adj_T_results.txt","adj_H_results.txt");
      //cout << adj_det.Nt <<endl;
      //cout << "displaying pointers" << endl;
      //cout << adj_det.steady_H << " " << adj_det.steady_T+11 << " " << src_adj.cumul_energies+100 <<endl;
    }

  double TEST = 0;
  if (strcmp(SIMU_TYPE.c_str(),"BOTH")==0 || strcmp(SIMU_TYPE.c_str(),"FORWARD")==0 || strcmp(SIMU_TYPE.c_str(),"REGULAR")==0 || strcmp(SIMU_TYPE.c_str(),"DIRECT")==0)
    {

      H_detect.show();
      cout << endl;
      T_detect.show();
      cout << endl << "STARTING FORWARD SIMULATION" << endl;
      // LOOP OVER PARTICLES
      for (int i=0; i<src.Npart; i++){
        print_percent(i, src.Npart);

	src.emit(&part,&r);
	if (0<part.pt0.x && part.pt0.x<1e-7 && 0<part.pt0.y && part.pt0.y<1e-7)
	  TEST = TEST + part.weight/Si.C/1e-14;
	//     cout << part.t << " " << src.t_max << endl;
	while (part.alive){
	  part.initiate_move(&r);

	  part.move(&ref, &presc, &per);

	  H_detect.measure(&part);
	  T_detect.measure(&part, &Si);

	  part.finish_move(&Si,&r, &ref, &presc, &per);

	}

      }
      // RECORD RESULTS
      H_detect.write("results_H.txt");
      T_detect.write("results_T.txt");
      // DISPLAYING RESULTS
      H_detect.show_results();
      T_detect.show_results();
    }
  return 0;


}
/*

{
  //INITIALIZE INPUTS
  time_t seconds;
  seconds = time(NULL);
  long seed = (long)seconds;
  
  RandomClass r;
  r.initialize(seed);
  
  int NPARTICLES, MAXIMUM_COLL;
  double TEQ;
  ifstream file;
  file.open("input.txt");
  if (file.fail()){
    cout << "NO INPUT FILE" << endl;
  }
  else
    {
      file >> NPARTICLES;
      file >> MAXIMUM_COLL;
      file >> TEQ;
    }
  // READ GEOMETRY FILES
  reflective_bdrs ref("reflective.txt");
  ref.show();
  prescribed_bdrs presc("prescribed.txt");
  presc.show();
  periodic_bdrs per("periodic.txt");
  per.show();

  // READ MATERIAL PROPERTIES
  materials Si("dataSi.txt", TEQ);
  Si.show_all();

  // READ SOURCE FILES
  sources src(&Si, &presc, "volumetric.txt", "body_force.txt", "initial.txt", NPARTICLES);
  src.display_type();
  cout << src.Np << endl;
  cout << src.Nv << endl;
  cout << src.Nb << endl;
  cout << src.Ni << endl;
  src.display_vol();
  src.display_bod();
  src.display_init();

  particle part(MAXIMUM_COLL);

  // DECLARE AND READ THE HEAT FLUX AND TEMPERATURE DETECTORS
  detector_array_H H_detect("H_detectors.txt", "times.txt");
  detector_array_T T_detect("T_detectors.txt", "times.txt");
  matlab_write_geometry(&ref, &presc, &per, &src, &T_detect, &H_detect, "geometry.m");
  H_detect.show();
  cout << endl;
  T_detect.show();

  // LOOP OVER PARTICLES
  for (int i=0; i<src.Npart; i++){
    print_percent(i, src.Npart);

    src.emit(&part,&r);
//    cout << part.t << " " << src.t_max << endl;
    while (part.alive){
      part.initiate_move(&r);
      part.move(&ref, &presc, &per);
      H_detect.measure(&part);
      T_detect.measure(&part, &Si);
      part.finish_move(&Si,&r, &ref, &presc, &per);

    }

  }
  // DISPLAY RESULTS IN TERMINAL AND RECORD THEM
  H_detect.show_results();
  H_detect.write("results_H.txt");
  T_detect.show_results();
  T_detect.write("results_T.txt");
  
  return 0;



}
*/
