// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef EDDA_READER
#define EDDA_READER

#include <string>
#include "edda_export.h"
#include "dataset/dataset.h"

namespace edda{
	
	/// \brief Read edda dataset from disk
	/// \param edda_file the name of the file to be read
	/// \return an object of Dataset, storing all info read from the file
	std::shared_ptr<Dataset<Real> > EDDA_EXPORT loadEddaScalarDataset(const std::string &edda_file);

	/// \brief Historical function to be deprecated. Instead, use loadEddaScalarDataset(const std::string &edda_file)
	/// \param edda_file the name of the file to be read
	/// \return an object of Dataset, storing all info read from the file
	std::shared_ptr<Dataset<Real> > EDDA_EXPORT loadEddaScalarDataset_noneVTK(const std::string &edda_file);


}

#endif
