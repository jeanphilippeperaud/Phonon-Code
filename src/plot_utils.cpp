

#include "reflective_bdrs.h"
#include "prescribed_bdrs.h"
#include "periodic_bdrs.h"
#include "sources.h"
#include "detector_array_H.h"
#include "detector_array_T.h"
#include <fstream>
#include "body_force.h"
#include "volumetric.h"
#include "initial.h"

void matlab_write_geometry(reflective_bdrs * ref, prescribed_bdrs * presc, periodic_bdrs * per, sources * src, detector_array_T * Td, detector_array_H * Hd, const char * filename)
{// this function writes a .m file that may be run in matlab to draw the geometry of the problem of interest
  // open the target file
  ofstream file;
  file.open(filename);
  file << "REFdata = zeros(" << ref->N << ",5);" << endl;
  file << "PREdata = zeros(" << presc->N << ",5);" << endl;
  file << "PERdata = zeros(" << per->N << ",6);" << endl;
  file << "BDdata = zeros(" << src->Nb << ",10);" << endl;
  file << "VOLdata = zeros(" << src->Nv << ",9);" << endl;
  file << "INITdata = zeros(" << src->Ni << ",9);" << endl;
  file << "T_detectors = zeros(" << Td->N << ",8);" << endl;
  file << "H_detectors = zeros(" << Hd->N << ",10);" << endl;

  // read the reflective boundaries
  for (int i=0; i<ref->N; i++)
    {
      file << "REFdata(" << i+1 << ",:)= [ " << ref->ref_handle[i].point1.x << " " << ref->ref_handle[i].point1.y  <<
	" " << ref->ref_handle[i].point2.x << " " << ref->ref_handle[i].point2.y  << " " << ref->ref_handle[i].spec << "];" << endl;
    }
  // read prescribed temperature boundaries
  for (int i=0; i<presc->N; i++)
    {
      file << "PREdata(" << i+1 << ",:)= [ " << presc->presc_handle[i].point1.x << " " << presc->presc_handle[i].point1.y  <<
	" " << presc->presc_handle[i].point2.x << " " << presc->presc_handle[i].point2.y  << " " << presc->presc_handle[i].temp-presc->presc_handle[i].temp_eq << "];" << endl;
    }
  // read periodic boundaries
  for (int i=0; i<per->N; i++)
    {
      file << "PERdata(" << i+1 << ",:)= [ " << per->per_handle[i].point1.x << " " << per->per_handle[i].point1.y  <<
	" " << per->per_handle[i].point2.x << " " << per->per_handle[i].point2.y  << " " << per->per_handle[i].vctr.x << " " << per->per_handle[i].vctr.y << "];" << endl;
    }
  // read source terms
  for (int i = 0; i<src->Nb; i++)
    {
      file << "BDdata(" << i+1 << ",:)= [ " << src->bd_frc_array[i].segs[0].point1.x << " " << src->bd_frc_array[i].segs[0].point1.y  << " " <<
	src->bd_frc_array[i].segs[1].point1.x << " " << src->bd_frc_array[i].segs[1].point1.y  << " " <<
	src->bd_frc_array[i].segs[2].point1.x << " " << src->bd_frc_array[i].segs[2].point1.y  << " " <<
	src->bd_frc_array[i].segs[3].point1.x << " " << src->bd_frc_array[i].segs[3].point1.y  << " " <<
	src->bd_frc_array[i].vctr.x << " " << src->bd_frc_array[i].vctr.y << "];" << endl;
    }
  for (int i = 0; i<src->Nv; i++)
    {
      file << "VOLdata(" << i+1 << ",:)= [ " << src->sprd_src_array[i].segs[0].point1.x << " " << src->sprd_src_array[i].segs[0].point1.y  << " " <<
	src->sprd_src_array[i].segs[1].point1.x << " " << src->sprd_src_array[i].segs[1].point1.y  << " " <<
	src->sprd_src_array[i].segs[2].point1.x << " " << src->sprd_src_array[i].segs[2].point1.y  << " " <<
	src->sprd_src_array[i].segs[3].point1.x << " " << src->sprd_src_array[i].segs[3].point1.y  << " " <<
	src->sprd_src_array[i].temp-src->sprd_src_array[i].temp_eq << "];" << endl;
    }
  for (int i = 0; i<src->Ni; i++)
    {
      file << "INITdata(" << i+1 << ",:)= [ " << src->initial_condition[i].segs[0].point1.x << " " << src->initial_condition[i].segs[0].point1.y  << " " <<
	src->initial_condition[i].segs[1].point1.x << " " << src->initial_condition[i].segs[1].point1.y  << " " <<
	src->initial_condition[i].segs[2].point1.x << " " << src->initial_condition[i].segs[2].point1.y  << " " <<
	src->initial_condition[i].segs[3].point1.x << " " << src->initial_condition[i].segs[3].point1.y  << " " <<
	src->initial_condition[i].temp-src->initial_condition[i].temp_eq << "];" << endl;
    }
  for (int i = 0; i<Td->N; i++)
    {
      file << "T_detectors(" << i+1 << ",:)= [ " << Td->t_handle[i].segs[0].point1.x << " " << Td->t_handle[i].segs[0].point1.y  << " " <<
	Td->t_handle[i].segs[1].point1.x << " " << Td->t_handle[i].segs[1].point1.y  << " " <<
	Td->t_handle[i].segs[2].point1.x << " " << Td->t_handle[i].segs[2].point1.y  << " " <<
	Td->t_handle[i].segs[3].point1.x << " " << Td->t_handle[i].segs[3].point1.y  << " " <<
	"];" << endl;
    }
  for (int i = 0; i<Hd->N; i++)
    {
      file << "H_detectors(" << i+1 << ",:)= [ " << Hd->h_handle[i].segs[0].point1.x << " " << Hd->h_handle[i].segs[0].point1.y  << " " <<
	Hd->h_handle[i].segs[1].point1.x << " " << Hd->h_handle[i].segs[1].point1.y  << " " <<
	Hd->h_handle[i].segs[2].point1.x << " " << Hd->h_handle[i].segs[2].point1.y  << " " <<
	Hd->h_handle[i].segs[3].point1.x << " " << Hd->h_handle[i].segs[3].point1.y  << " " <<
	Hd->h_handle[i].vctr.x << " " << Hd->h_handle[i].vctr.y  << " " << "];" << endl;
    }
  //code for drawing
  file << "figure; hold on;" << endl;
  file << "for i=1:" << src->Nb << endl;
  file << "hb(i) = fill(BDdata(i,1:2:8),BDdata(i,2:2:8),'g','EdgeColor','none');" << endl;
  file << "end" << endl;
  file << "for i=1:" << src->Nv << endl;
  file << "hv(i) = fill(VOLdata(i,1:2:8),VOLdata(i,2:2:8),'r',LineWidth',2);" << endl;
  file << "end" << endl;
  file << "for i=1:" << src->Ni << endl;
  file << "hi(i) = fill(INITdata(i,1:2:8),INITdata(i,2:2:8),'y','EdgeColor','none');" << endl;
  file << "end" << endl;
  file << "for i=1:" << Td->N << endl;
  file << "ht(i) = fill(T_detectors(i,1:2:8),T_detectors(i,2:2:8),'w','FaceColor','none');" << endl;
  file << "set(ht(i),'LineStyle','--');" << endl;
  file << "end" << endl;
  file << "for i=1:" << Hd->N << endl;
  file << "hh(i) = fill(H_detectors(i,1:2:8),H_detectors(i,2:2:8),'w','FaceColor','none');" << endl;
  file << "set(hh(i),'LineStyle','--');" << endl;
  file << "end" << endl;

  file << "for i=1:" << ref->N << endl;
  file << "plot(REFdata(i,1:2:3),REFdata(i,2:2:4),'k','LineWidth',2);" << endl;
  file << "end" << endl;
  file << "for i=1:" << presc->N << endl;
  file << "plot(PREdata(i,1:2:3),PREdata(i,2:2:4),'b','LineWidth',2);" << endl;
  file << "end" << endl;
  file << "for i=1:" << per->N << endl;
  file << "plot(PERdata(i,1:2:3),PERdata(i,2:2:4),'--','LineWidth',2);" << endl;
  file << "end" << endl;
}
