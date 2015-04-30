// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef COMMON_H_
#define COMMON_H_

#include <cstdint>

#include "edda.h"

namespace edda {

const double DEG_TO_RAD = 0.0174532925199432957692;  // PI / 180
const double RAD_TO_DEG = 57.2957795130823208768;  // 180 / PI
const double PIBY2 = 1.57079632679489661923;  // PI / 2
const double EPS = 1.0e-6;

enum ReturnStatus { SUCCESS = 0, FAIL };

#ifdef OS_WIN
typedef long long int64_t;
#endif

// May not be needed in the release version but just keep here for now.
class NotImplementedException {};

}  // namespace edda

#endif  // COMMON_H_
