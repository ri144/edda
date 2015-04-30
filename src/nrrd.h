#ifndef NRRD_H
#define NRRD_H

#include "common.h"

namespace edda{
ReturnStatus write_nrrd_3d(const char *nrrd_fname, const char *raw_fname, int w, int h, int d, const char *type);

}

#endif

