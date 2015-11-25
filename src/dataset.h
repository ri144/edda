// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DATASET_H_
#define DATASET_H_

#include <vector>
#include <iostream>
#include <memory>
#include "geometry/grid.h"
#include "core/data_array.h"
#include <distributions/distribution.h>

namespace edda {

template <typename T>  // Return type of at_phys
class Dataset {
protected:
    Grid *pGrid;
    AbstractDataArray *pArray;
public:
    Dataset(Grid *pGrid, AbstractDataArray *pArray) {
        this->pGrid = pGrid;
        this->pArray = pArray;
    }
    ~Dataset() {
        delete pGrid;
        delete pArray;
    }

    Grid *getGrid() { return pGrid; }
    AbstractDataArray *getArray () {return pArray; }

    int *getDimension() {
      CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(pGrid) ;
      if (!cartesianGrid)
        throw NotImplementedException();
      return cartesianGrid->getDimension();
    }

    ReturnStatus at_phys(const VECTOR3 &pos, T &output) const {
        ReturnStatus r;
        switch (pGrid->getInterpType())
        {
        case InterpType::TRI_LERP:
        {
            PointInfo pinfo;
            pinfo.phyCoord = pos ;
            r = pGrid->phys_to_cell(pinfo);
            if (r != SUCCESS)
                return r;

            std::vector<int> vVertices;
            r = pGrid->getCellVertices(pinfo.inCell, vVertices);
            if (r != SUCCESS)
                return r;

            int i;
            std::vector<T> vData(vVertices.size());
            for (i=0; i<(int)vVertices.size(); i++)
                vData[i] = boost::any_cast<T> ( pArray->getItem(vVertices[i]) );

            output = triLerp(vData[0], vData[1], vData[2], vData[3],
                             vData[4], vData[5], vData[6], vData[7],
                             pinfo.interpolant.getData());
            break;
        }
        default:
            throw NotImplementedException();
            break;
        }

        return SUCCESS;
    }

    // value at computational space, only for structured grids
    const T &at_comp(int i, int j, int k) const {
      ReturnStatus r;
      CartesianGrid *cartesianGrid = dynamic_cast<CartesianGrid *>(pGrid) ;

      if (cartesianGrid) {
        int idx;
        r = cartesianGrid->getIndex(i,j,k, idx);
        if (r!=SUCCESS)
          throw OutOfBoundException();

        return boost::any_cast<const T&>( pArray->getItem(idx) );
      }

      // TODO for other grid types
      throw NotImplementedException();
    }

};


// create a shared pointer of Dataset
template <typename T>  // Return type of at_phys
inline std::shared_ptr< Dataset<T> >
make_Dataset(Grid *pGrid, AbstractDataArray *pArray)
{
    return std::shared_ptr< Dataset<T> >(
                new Dataset<T> (pGrid, pArray) );
}

} // namespace edda

#endif // STRUCTURED_GRID_DATASET_H_
