/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vector Field
//                 Created: Han-Wei Shen, Liya Li
//                 The Ohio State University
//                 Date:		06/2005
//                 Vector Field: 3D Static or Time-Varying
//
///////////////////////////////////////////////////////////////////////////////

// Field wraps solutions and grids

#ifndef _VECTOR_FIELD_H_
#define _VECTOR_FIELD_H_

#include <assert.h>
#include "header.h"
#include "VectorMatrix.h"
#include "Grid.h"
#include "Solution.h"

namespace edda{

//////////////////////////////////////////////////////////////////////////
// field class
//////////////////////////////////////////////////////////////////////////
template<class PointType>
class Field
{
public:
    // constructor and destructor
    Field()  {}
    virtual ~Field() {}

    //////////////////////////////////////////////////////////////////////////
    // to obtain one or more field values at cell vertices
    //
    // input
    //		cellId:			from which cell to obtain value
    //		t:				which time step
    // output
    //		vNodeData:	include one or more node values
    //////////////////////////////////////////////////////////////////////////
    virtual ReturnStatus at_cell(int cellId, const double t, std::vector<VECTOR3>& vNodeData)=0;

    //////////////////////////////////////////////////////////////////////////
    // to obtain node data at the physical position pos in timestep t
    //
    // input
    //		pos:		physical position of node
    //		t:			time step
    //		fromCell:	if not -1, which cell this position is generated from
    // output
    //		nodeData:	data value of this node
    //      pInfo:      point information in the cell
    //////////////////////////////////////////////////////////////////////////
    virtual ReturnStatus at_phys(const int fromCell, VECTOR3& pos, PointInfo& pInfo,const float t, PointType& nodeData)=0;
    virtual ReturnStatus at_phys(const VECTOR3 &pos, float t, PointType& nodeData)=0;
    virtual float volume_of_cell(int cellId)=0;

    virtual void scan (
        void (*func)(int iLocalT, int iNode, VECTOR3 *pv3)
       )=0;

    virtual int getTimeSteps(void) =0;
    virtual int getMinTimeStep(void) =0;
    virtual int getMaxTimeStep(void) =0;
    virtual void getBoundary(VECTOR3& minB, VECTOR3& maxB) =0;
    virtual void setBoundary(VECTOR3 minB, VECTOR3 maxB) =0;
    virtual bool isInRealBoundaries(PointInfo& p, float time=0) =0;

};

//////////////////////////////////////////////////////////////////////////
// field class
//////////////////////////////////////////////////////////////////////////
template<class PointType>
class GeneralField: public Field<PointType>
{
private:
	Grid* m_pGrid;						// grid
    Solution<PointType>* m_pSolution;				// vector data
	int m_nTimeSteps;
	bool m_bIsNormalized;				// whether the solution is normalized or not
	int m_MinT, m_MaxT; // the min and max time step of the data field

public:
	// constructor and destructor
    GeneralField() :
        m_nTimeSteps(0),
        m_pGrid (NULL),
        m_pSolution (NULL),
        m_bIsNormalized (false),
        m_MinT(-1),
        m_MaxT(-1)
    {}

    GeneralField(Grid* pGrid, Solution<PointType>* pSolution, int timesteps, int min_t=0) :
        m_pGrid (pGrid),
        m_pSolution (pSolution),
        m_nTimeSteps (timesteps),
        m_bIsNormalized (false),
        m_MinT (min_t),
        m_MaxT (min_t + timesteps -1)
    {
        assert((pGrid != NULL) && (pSolution != NULL));
    }

    virtual ~GeneralField() {
        if(m_pGrid != NULL)	delete m_pGrid;
        if(m_pSolution != NULL) delete m_pSolution;
    }

    inline CellType getCellType(void) { return m_pGrid->GetCellType(); }


    virtual ReturnStatus at_cell(int cellId, const double t, std::vector<VECTOR3>& vNodeData)
    {
        VECTOR3 nodeData;
        std::vector<int> vVerIds;
        int iFor;

        if(m_pGrid->getCellVertices(cellId, vVerIds) == FAIL)
            return FAIL;

        for(iFor = 0; iFor < (int)vVerIds.size(); iFor++)
        {
            if(m_pSolution->GetValue(vVerIds[iFor], t, nodeData) == FAIL)
                return FAIL;
            vNodeData.push_back(nodeData);
        }
        return SUCCESS;
    }

    virtual ReturnStatus at_phys(const VECTOR3 &pos, float t, PointType& data)
    {
        PointInfo pInfo;
        return at_phys(-1, pos, pInfo, t, data);
    }


    virtual int at_phys(const int fromCell, VECTOR3& pos, PointInfo& pInfo,const float t, PointType& data)
    {
        std::vector<PointType> vNodeData;

        // find the cell this position belongs to
        pInfo.phyCoord = pos;
        pInfo.fromCell = fromCell;
        pInfo.inCell = -1;
        if(m_pGrid->phys_to_cell(pInfo) == FAIL) {
            return FAIL;
        }

        // get vertex value at cell vertices
        if(at_cell(pInfo.inCell, t, vNodeData) == FAIL) {
            return FAIL;
        }

        // interpolate in the cell
        m_pGrid->interpolate(data, vNodeData, pInfo.interpolant);

        return SUCCESS;
    }

    //virtual int at_comp(const int i, const int j, const int k, const float t, VECTOR3& dataValue);
    virtual float volume_of_cell(int cellId) {
        return m_pGrid->cellVolume(cellId);
    }

	virtual void Scan
	  (
	   void (*func)(int iLocalT, int iNode, VECTOR3 *pv3)
      )
    {
        m_pSolution->Scan(func);
    }

    virtual int getTimeSteps(void) {return m_nTimeSteps;}
    virtual int getMinTimeStep(void) {return m_MinT;}
    virtual int getMaxTimeStep(void) {return m_MaxT;}
    virtual void getBoundary(VECTOR3& minB, VECTOR3& maxB) { m_pGrid->Boundary(minB, maxB); }
    virtual void setBoundary(VECTOR3 minB, VECTOR3 maxB) {m_pGrid->SetBoundary(minB, maxB);	}
    virtual void BoundaryIntersection(VECTOR3& intersectP,VECTOR3& startP,VECTOR3& endP,float* stepSize,float oldStepSize)
        { m_pGrid->BoundaryIntersection(intersectP, startP, endP, stepSize, oldStepSize); }
    virtual bool isInRealBoundaries(PointInfo& p, float time=0)
    {
        return m_pGrid->isInRealBBox(p.phyCoord, time);
    }

    // specific functions

    inline Grid *getGrid() {return m_pGrid;}
    inline Solution<PointType> *getSolution() {return m_pSolution; }

protected:
	// field functions
    virtual bool isTimeVarying(void)
    {
        return (m_nTimeSteps > 1);
    }

};

template<class PointType>
class StructuredField: public GeneralField<PointType>
{
public:
    // specific functions
    //////////////////////////////////////////////////////////////////////////
    // to obtain node data at position (i, j, k) in time t
    //
    // input
    //		(i, j, k):	position in computational space
    //		t:			time step
    // output
    //		dataValue:	data value of this position
    // return
    //		1:			operation successful
    //		-1:			operation fail
    //////////////////////////////////////////////////////////////////////////
    virtual ReturnStatus at_vert(const int i, const int j, const int k, const float t, VECTOR3& dataValue)
    {
        int xdim, ydim, zdim;
        getDimension(xdim, ydim, zdim);
        if (i<0 || j<0 || k<0 || i>=xdim || j>=ydim || k>=zdim)
            return -1;
        return this->m_pSolution->GetValue((k*ydim*xdim+j*xdim+i), t, dataValue);
    }

    // get the dimension in computational space (structured grids only)
    virtual void getDimension(int& xdim, int& ydim, int& zdim);
};


template<class PointType>
class UnstructuredField: public GeneralField<PointType>
{
public:
};

} // namespace edda

#endif
