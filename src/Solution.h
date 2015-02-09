/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vector Field
//                 Created: Han-Wei Shen, Liya Li
//                 The Ohio State University
//                 Date:		06/2005
//                 Vector Field Data
//
///////////////////////////////////////////////////////////////////////////////



#ifndef _SOLUTION_H_
#define _SOLUTION_H_

#include <cmath>
#include <cstdlib>
#include <cfloat>
#include <cassert>
#include "header.h"
#include "VectorMatrix.h"
#include "Interpolator.h"

namespace edda{

template<class PointType>
class Solution
{
private:
    PointType** m_pDataArray;				// data value
    PointType* m_pMinValue;				// value with min magnitude for each time step
    PointType* m_pMaxValue;				// value with max magnitude for each time step
	int m_nNodeNum;						// how many nodes each time step
	int m_nTimeSteps;					// how many time steps
	float m_fMinMag;					// minimal magnitude
	float m_fMaxMag;					// maximum magnitude
	int m_MinT, m_MaxT;

public:
	// constructor
    Solution() {Reset();}
    //////////////////////////////////////////////////////////////////////////
    // input
    //		pData:		this is a 2D array for storing static or
    //					time-varying data
    //		nodeNum:	number of total nodes
    //		timeSteps:	how many time steps, for static data, this is 1
    //////////////////////////////////////////////////////////////////////////
    Solution(PointType** pData, int nodeNum, int timeSteps){
        assert((pData != NULL) && (nodeNum > 0) && (timeSteps > 0));

        m_nNodeNum = nodeNum;
        m_nTimeSteps = timeSteps;
        m_MinT = 0; m_MaxT = timeSteps-1;
        m_pDataArray = pData;
    }
    Solution(PointType** pData, int nodeNum, int timeSteps, int min_t, int max_t){
        assert((pData != NULL) && (nodeNum > 0) && (timeSteps > 0));

        m_nNodeNum = nodeNum;
        m_nTimeSteps = timeSteps;
        m_MinT = min_t; m_MaxT = max_t;
        m_pDataArray = pData;
    }

    ~Solution(){
        int iFor;

        if(m_pDataArray != NULL)
        {
            for(iFor = 0; iFor < m_nTimeSteps; iFor++)
            {
                delete[] m_pDataArray[iFor];
            }
            delete[] m_pDataArray;
        }
    }

	void Reset();

	// solution functions
	void SetMinMaxTime(int min_t, int max_t) {m_MinT = min_t; m_MaxT = max_t;}
    void SetValue(int t, PointType* pData, int nodeNum);
    int GetMinMaxValueAll(PointType& minVal, PointType& maxVal);
    int GetMinMaxValue(int t, PointType& minVal, PointType& maxVal);
	void ComputeMinMaxValue(void);
	bool isTimeVarying(void);
    ReturnStatus GetValue(int id, const float t, PointType& nodeData);
	void Normalize(bool bLocal);
	void Scale(float scaleF);
    void Translate(PointType& translate);
	// ADD-BY-LEETEN 02/02/2012-BEGIN
	// This function scan the solution with the user-specified function func().
	void
	  Scan
	  (
	   void (*func)(int iLocalT, int iNode, VECTOR3* pv3)
	   );
	// ADD-BY-LEETEN 02/02/2012-END
};



template<class PointType>
void Solution<PointType>::Reset()
{
    m_pDataArray = NULL;
    m_nNodeNum = 0;
    m_nTimeSteps = 1;
    m_MinT = 0; m_MaxT = 0;
}

//////////////////////////////////////////////////////////////////////////
// change vector data on-the-fly
//////////////////////////////////////////////////////////////////////////
template<class PointType>
void Solution<PointType>::SetValue(int t, PointType* pData, int nodeNum)
{
    assert( t>=m_MinT && t<=m_MaxT );

    t = t-m_MinT;
    m_pDataArray[t] = new VECTOR3[nodeNum];
    assert(m_pDataArray[t] != NULL);
    for(int jFor = 0; jFor < nodeNum; jFor++)
        m_pDataArray[t][jFor] = pData[jFor];

}

//////////////////////////////////////////////////////////////////////////
// whether field is time varying
//////////////////////////////////////////////////////////////////////////
template<class PointType>
bool Solution<PointType>::isTimeVarying(void)
{
    return (m_nTimeSteps > 1);
}

//////////////////////////////////////////////////////////////////////////
// get value of node id at time step t
// input
//		id:			node Id
//		t:			time step in check
// output
//		nodeData:	vector value at this node
// return
//		1:			operation successful
//		-1:			invalid id
//////////////////////////////////////////////////////////////////////////
template<class PointType>
ReturnStatus Solution<PointType>::GetValue(int id, float t, PointType& nodeData)
{
  float adjusted_t = t - m_MinT;
    if((id < 0) || (id >= m_nNodeNum) || (adjusted_t < 0.0) || (adjusted_t > (float)(m_nTimeSteps-1)))
        return FAIL;

    if(!isTimeVarying())
        nodeData = m_pDataArray[(int)adjusted_t][id];
    else
    {
        int lowT, highT;
        float ratio;
        lowT = (int)floor(adjusted_t);
        ratio = adjusted_t - (float)floor(adjusted_t);
        highT = lowT + 1;
        if(lowT >= (m_nTimeSteps-1))
        {
            highT = lowT;
            ratio = 0.0;
        }
        nodeData.Set(Lerp(m_pDataArray[lowT][id][0], m_pDataArray[highT][id][0], ratio),
                     Lerp(m_pDataArray[lowT][id][1], m_pDataArray[highT][id][1], ratio),
                     Lerp(m_pDataArray[lowT][id][2], m_pDataArray[highT][id][2], ratio));
    }

    return SUCCESS;
}

//////////////////////////////////////////////////////////////////////////
// to normalize the vector field
// input
// bLocal: whether to normalize in each timestep or through all timesteps.
//		   if locally, then divide its magnitude; if globally, then divide
//		   by the maximal magnitude through the whole field
//////////////////////////////////////////////////////////////////////////
template<class PointType>
void Solution<PointType>::Normalize(bool bLocal)
{
    int iFor, jFor;
    float mag, u, v, w;

    m_fMinMag = FLT_MAX;
    m_fMaxMag = -FLT_MAX;

    if(bLocal)
    {
        for(iFor = 0; iFor < m_nTimeSteps; iFor++)
        {
            for(jFor = 0; jFor < m_nNodeNum; jFor++)
            {
                mag = m_pDataArray[iFor][jFor].GetMag();
                if(mag != 0.0)
                {
                    u = m_pDataArray[iFor][jFor][0]/mag;
                    v = m_pDataArray[iFor][jFor][1]/mag;
                    w = m_pDataArray[iFor][jFor][2]/mag;
                    m_pDataArray[iFor][jFor].Set(u, v, w);
                }

                if(mag < m_fMinMag)
                    m_fMinMag = mag;
                if(mag > m_fMaxMag)
                    m_fMaxMag = mag;
            }
        }
    }
    else
    {
        for(iFor = 0; iFor < m_nTimeSteps; iFor++)
        {
            for(jFor = 0; jFor < m_nNodeNum; jFor++)
            {
                mag = m_pDataArray[iFor][jFor].GetMag();
                if(mag < m_fMinMag)
                    m_fMinMag = mag;
                if(mag > m_fMaxMag)
                    m_fMaxMag = mag;
            }
        }

        for(iFor = 0; iFor < m_nTimeSteps; iFor++)
        {
            for(jFor = 0; jFor < m_nNodeNum; jFor++)
            {
                u = m_pDataArray[iFor][jFor][0]/m_fMaxMag;
                v = m_pDataArray[iFor][jFor][1]/m_fMaxMag;
                w = m_pDataArray[iFor][jFor][2]/m_fMaxMag;
                m_pDataArray[iFor][jFor].Set(u, v, w);
            }
        }
    }
}

template<class PointType>
void Solution<PointType>::Scale(float scaleF)
{
    int iFor, jFor;
    float u, v, w;

    for(iFor = 0; iFor < m_nTimeSteps; iFor++)
      {
        for(jFor = 0; jFor < m_nNodeNum; jFor++)
          {
        u = m_pDataArray[iFor][jFor][0]*scaleF;
        v = m_pDataArray[iFor][jFor][1]*scaleF;
        w = m_pDataArray[iFor][jFor][2]*scaleF;
        m_pDataArray[iFor][jFor].Set(u, v, w);
          }
      }
}

template<class PointType>
void Solution<PointType>::Translate(PointType& translate) {
    int iFor, jFor;
    float u, v, w;

    for (iFor = 0; iFor < m_nTimeSteps; iFor++) {
        for (jFor = 0; jFor < m_nNodeNum; jFor++) {
            u = m_pDataArray[iFor][jFor][0] + translate[0];
            v = m_pDataArray[iFor][jFor][1] + translate[1];
            w = m_pDataArray[iFor][jFor][2] + translate[2];
            m_pDataArray[iFor][jFor].Set(u, v, w);
        }
    }
}

// ADD-BY-LEETEN 02/02/2012-BEGIN
template<class PointType>
void Solution<PointType>::Scan
(
 void (*_func)(int iLocalT, int iNode, VECTOR3* pv3)
)
{
    int iFor, jFor;
    float u, v, w;

    for(int iFor = 0; iFor < m_nTimeSteps; iFor++)
      {
        for(int jFor = 0; jFor < m_nNodeNum; jFor++)
          {
        _func
          (
              iFor,
              jFor,
              &m_pDataArray[iFor][jFor]
          );
          }
      }
}
// ADD-BY-LEETEN 02/02/2012-END

// compute the min and max value with minimal and maximal magnitude
template<class PointType>
void Solution<PointType>::ComputeMinMaxValue(void)
{
    int indexMin, indexMax, iFor, jFor;
    float minMag, maxMag, mag;

    m_pMinValue = new VECTOR3[m_nTimeSteps];
    m_pMaxValue = new VECTOR3[m_nTimeSteps];

    for(iFor = 0; iFor < m_nTimeSteps; iFor++)
    {
        minMag = FLT_MAX;
        maxMag = -FLT_MAX;
        for(jFor = 0; jFor < m_nNodeNum; jFor++)
        {
            mag = m_pDataArray[iFor][jFor].GetMag();
            if(mag < minMag)
            {
                minMag = mag;
                indexMin = jFor;
            }
            if(mag > maxMag)
            {
                maxMag = mag;
                indexMax = jFor;
            }
        }

        m_pMinValue[iFor] = m_pDataArray[iFor][indexMin];
        m_pMaxValue[iFor] = m_pDataArray[iFor][indexMax];
    }
}

// get the min and max value for all time steps
template<class PointType>
int Solution<PointType>::GetMinMaxValueAll(PointType& minVal, PointType& maxVal)
{
    minVal = m_pMinValue[0];
    maxVal = m_pMaxValue[0];

    for(int tFor = 1; tFor < m_nTimeSteps; tFor++)
    {
        if(minVal.GetMag() > m_pMinValue[tFor].GetMag())
            minVal = m_pMinValue[tFor];

        if(maxVal.GetMag() < m_pMaxValue[tFor].GetMag())
            maxVal = m_pMaxValue[tFor];
    }

    return 1;
}

// get the min and max value for timestep t
template<class PointType>
int Solution<PointType>::GetMinMaxValue(int t, PointType& minVal, PointType& maxVal)
{
        t = t - m_MinT;
    if((t >= 0) && (t < m_nTimeSteps))
    {
        minVal = m_pMinValue[t];
        maxVal = m_pMaxValue[t];
        return 1;
    }
    else
        return -1;
}



} // namespace edda

#endif
