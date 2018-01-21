// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef EDDA_WRITER
#define EDDA_WRITER

#include <string>
#include "edda_export.h"
#include "dataset/dataset.h"

namespace edda{

	/// \brief Save an object of edda Dataset to disk
	/// \param dataset the object of edda Dataset to be saved
	/// \param edda_file the name of the file to be saved to
	void EDDA_EXPORT writeEddaDataset(std::shared_ptr<Dataset<Real> > dataset, const std::string &edda_file);
}

#endif
