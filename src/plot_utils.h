


#ifndef PLOT_UTILS_H
#define PLOT_UTILS_H

class reflective_bdrs;
class periodic_bdrs;
class prescribed_bdrs;
class sources;
class detector_array_T;
class detector_array_H;

void matlab_write_geometry(reflective_bdrs * ref, prescribed_bdrs * presc, periodic_bdrs * per, sources * src, detector_array_T * Td, detector_array_H * Hd, const char * filename);

#endif
