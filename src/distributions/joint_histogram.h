// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DIST_JOINT_HISTOGRAM_H_
#define DIST_JOINT_HISTOGRAM_H_

#include <cstdlib>
#include <iostream>
#include <vector>
#include <unordered_set>
#define _USE_MATH_DEFINES  // For Visual Studio
#include <math.h>

#include <boost/limits.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>
#include <boost/bind.hpp>

#include "common.h"
#include "distribution_tag.h"
#include "core/statistics.h"
//#include "invert_matrix.h"
#include <boost/unordered_map.hpp>

namespace edda {
namespace dist {

// ------------------------------------------------------------------------------
///
/// \brief Defines a Joint Histogram class
///
struct EDDA_EXPORT JointHistogram: public DiscreteDistributionTag, public JointDistributionTag {
	ublas_vector mean; // boost's vector
	
	// constructor
	__host__ __device__
	JointHistogram() {
		// initialize a default joint histogram for 3 vars
		num_comps = 3;
	
		min_vals = std::vector<Real>(num_comps, 0);
		max_vals = std::vector<Real>(num_comps, 0);
		num_bins = std::vector<int>(num_comps, 0);
		mean = ublas::zero_vector<Real>(3);
		setCovariance( ublas::identity_matrix<Real>(3) );
	}

	__host__ __device__
	JointHistogram(std::vector<Real*>& dataAry, int nElements, const std::vector<Real> &mins, 
					const std::vector<Real> &maxs, const std::vector<int> &nBins) {
		assert(dataAry.size()==mins.size() && mins.size()==maxs.size() && mins.size()==nBins.size());
		for(int i=0; i<mins.size(); i++) assert(maxs[i]>=mins[i]);

		// initialize a joint histogram from given data
		num_comps = dataAry.size();
		min_vals = mins;
		max_vals = maxs;
		num_bins = nBins;
	
		// compute the pdf, mean and cov 
		_update(dataAry, nElements);
	}

	__host__ __device__
	JointHistogram(int comps, std::vector<Real> mins, std::vector<Real> maxs, std::vector<Real> binWidths, 
					std::vector<int> nBins, boost::unordered_map<std::vector<int>, Real> new_pdf, 
					ublas_vector new_mean, ublas_matrix new_cov) {
		// construct from a given JointHistogram parameters
		num_comps = comps;
		min_vals = mins;
		max_vals = maxs;
		num_bins = nBins;
		bin_widths = binWidths;
		pdf = new_pdf;
		mean = new_mean;
		setCovariance(new_cov);
		_computeJointCdf();
	}

	__host__ __device__
	void setMinVals(const std::vector<Real>& mins) { min_vals = mins; } 
	__host__ __device__
	std::vector<Real> getMinVals() const { return  min_vals; } 

	__host__ __device__
	void setMaxVals(const std::vector<Real>& maxs) { max_vals = maxs; } 
	__host__ __device__
	std::vector<Real> getMaxVals() const { return  max_vals; } 

	__host__ __device__
	void setBinWidths(const std::vector<Real>& widths) { bin_widths = widths; } 
	__host__ __device__
	std::vector<Real> getBinWidths() const { return  bin_widths; } 

	__host__ __device__
	void setNumBins(const std::vector<int>& bins) { num_bins = bins; } 
	__host__ __device__
	std::vector<int> getNumBins() const { return  num_bins; } 

	__host__ __device__
	void setNumComps(const int comps) { num_comps = comps; } 
	__host__ __device__
	int getNumComps() const { return  num_comps; } 

	boost::unordered_map<std::vector<int>, Real> getDistr() const { return pdf; }
	void setDistr(boost::unordered_map<std::vector<int>, Real> p) { pdf = p; }

	std::vector<std::pair<std::vector<int>, Real>> getJointCdf() const { 
		assert(cdf.size()>0); 
		assert(cdf.size()==pdf.size()); 
		return cdf; 
	}
	
	__host__ __device__
	void setCovariance(const ublas_matrix &cov) {
		this->cov = cov;
		bool r = invert_matrix(this->cov, inv_cov);
		if (!r)
			std::cout << "Error: inverse covariance cannot be computed for matrix: " << cov << std::endl;
		// for debugging
		//std::cout << "(JointHistogram debug) Inverse func: " << inv_cov << std::endl;
		//std::cout << "(JointHistogram debug) " << ublas::prod(cov, inv_cov) << std::endl;
		this->det = determinant(this->cov);
		//std::cout << "(JointHistogram debug) Det: " << det << std::endl;
	}
	
	__host__ __device__
	const ublas_matrix &getCovariance() const {  return this->cov ; }

	__host__ __device__
	const ublas_matrix &getInvCovariance() const {  return this->inv_cov ; }

	__host__ __device__
	double getDet() const {  return this->det ; }

	///
	///\brief project to which dimensions, i.e. vars
	///
	JointHistogram marginalization(std::unordered_set<int>& vars) const {
		assert(pdf.size()>0);
		assert(vars.size()>0);

		std::vector<int> new_num_bins;
		std::vector<Real> new_bin_widths;
		std::vector<Real> new_min_vals;
		std::vector<Real> new_max_vals;
		boost::unordered_map<std::vector<int>, Real> new_pdf;
		ublas_vector new_mean = ublas::zero_vector<Real>(vars.size()); // boost's vector
		ublas_matrix new_cov(ublas::identity_matrix<Real>(vars.size())); // boost's vector
		int varid = 0;
		for(auto itv=vars.begin(); itv!=vars.end(); ++itv) {
			new_num_bins.push_back(num_bins[*itv]);
			new_bin_widths.push_back(bin_widths[*itv]);
			new_min_vals.push_back(min_vals[*itv]);
			new_max_vals.push_back(max_vals[*itv]);
			new_mean[varid] = mean[*itv];
			varid++;
		}
		// marginalize the pdf
		for (auto itp=pdf.begin(); itp!=pdf.end(); ++itp) {
			std::vector<int> bin = itp->first;
			std::vector<int> key;
			for(auto itv=vars.begin(); itv!=vars.end(); ++itv) {
				key.push_back(bin[*itv]);
			}
			new_pdf[key] += itp->second;
		}
		// compute new mean and cov
		int iitr=0;	
		for(auto itv=vars.begin(); itv!=vars.end(); ++itv) {
			int jitr = 0;
			for(auto jtv=vars.begin(); jtv!=vars.end(); ++jtv) {
				new_cov(iitr,jitr) = cov(*itv, *jtv);
				jitr++;
			}
			iitr++;
		}
		JointHistogram jhist(vars.size(), new_min_vals, new_max_vals, new_bin_widths, new_num_bins, new_pdf, new_mean, new_cov);
		return jhist; 
	}
	
	JointHistogram conditionalHist(std::unordered_set<int>& vars, std::vector<int>& cond_var, std::vector<std::pair<int,int>>& bin_range) { 
		// step 0: check input parameter
		assert(pdf.size()>0);
		assert(vars.size()>0);
		assert(cond_var.size()>0 && cond_var.size()==bin_range.size());
		for(int i=0; i<bin_range.size(); i++) {// validate bin range
			assert(bin_range[i].first<bin_range[i].second);
		}
		for(int i=0; i<cond_var.size(); i++) {// check for duplication
			for(int j=i+1; j<cond_var.size(); j++) {
				assert(cond_var[i]!=cond_var[j]);
			}
		}
		// setp 1: computing P(cond_var in bin_range)
		Real pB = 0;
		bool flag_valid = false;
		for (auto itp=pdf.begin(); itp!=pdf.end(); ++itp) {
			std::vector<int> bin = itp->first;
			bool flag_in_range = true;
			for(int i=0; i<cond_var.size(); i++) {
				if(bin[cond_var[i]]<bin_range[i].first || bin[cond_var[i]]>bin_range[i].second) {
					flag_in_range = false;
					break;
				}
			}
			if(flag_in_range){
				flag_valid = true;
				pB += itp->second;
			}
		}
		if(pB==0 || !flag_valid) {// no point in conditional range
			throw std::runtime_error("no data item falls in the given bin ranges!");
		}
	
		// step 2: computing P(vars AND cond_vars)/P(cond_var in bin_range) 
		boost::unordered_map<std::vector<int>, Real> new_pdf;
		for (auto itp=pdf.begin(); itp!=pdf.end(); ++itp) {
			std::vector<int> bin = itp->first;
			std::vector<int> key;
			for(auto itv=vars.begin(); itv!=vars.end(); ++itv) {
				key.push_back(bin[*itv]);
			}
			bool flag_in_range = true;
			for(int i=0; i<cond_var.size(); i++) {
				if(bin[cond_var[i]]<bin_range[i].first || bin[cond_var[i]]>bin_range[i].second) {
					flag_in_range = false;
					break;
				}
			}
			if(flag_in_range)
				new_pdf[key] += itp->second/pB;
		}
		
		// step 3: construct a new JointHistogram
		std::vector<int> new_num_bins;
		std::vector<Real> new_bin_widths;
		std::vector<Real> new_min_vals;
		std::vector<Real> new_max_vals;
		ublas_vector new_mean = ublas::zero_vector<Real>(vars.size()); // boost's vector
		ublas_matrix new_cov(ublas::identity_matrix<Real>(vars.size())); // boost's vector
		int varid = 0;
		for(auto itv=vars.begin(); itv!=vars.end(); ++itv) {
			new_num_bins.push_back(num_bins[*itv]);
			new_bin_widths.push_back(bin_widths[*itv]);
			new_min_vals.push_back(min_vals[*itv]);
			new_max_vals.push_back(max_vals[*itv]);
			new_mean[varid] = mean[*itv];
			varid++;
		}
		int iitr=0;	
		for(auto itv=vars.begin(); itv!=vars.end(); ++itv) {
			int jitr = 0;
			for(auto jtv=vars.begin(); jtv!=vars.end(); ++jtv) {
				new_cov(iitr,jitr) = cov(*itv, *jtv);
				jitr++;
			}
			iitr++;
		}
		JointHistogram jhist(vars.size(), new_min_vals, new_max_vals, new_bin_widths, new_num_bins, new_pdf, new_mean, new_cov);
		return jhist; 
	}
	
private:
	// histogram information
	int num_comps;
	std::vector<Real> min_vals;
	std::vector<Real> max_vals;
	std::vector<int> num_bins;
	std::vector<Real> bin_widths;
	boost::unordered_map<std::vector<int>, Real> pdf;
	std::vector<std::pair<std::vector<int>, Real>> cdf;
	// cov matrix 
	ublas_matrix cov;
	ublas_matrix inv_cov;
	double det;

	// compute joint PDF, mean, and cov 
	__host__ __device__
	void _update(std::vector<Real*>& dataAry, int nElements) {
		// 1. compute the PDF (i.e. a sparce matrix)
		boost::unordered_map<std::vector<int>, int> cntMap;
		std::vector<int> key(num_comps);
		bin_widths.clear();	
		for(int i=0; i<num_comps; i++)	
			bin_widths.push_back((max_vals[i]-min_vals[i])/num_bins[i]);
		for(int j=0; j<nElements; j++) {
			for(int i=0; i<num_comps; i++) {// value to bin index
				key[i] = std::floor((dataAry[i][j]-min_vals[i])/bin_widths[i]);
				if(key[i]>=0 && key[i]<num_bins[i])
					continue;	
				else
					key[i] = num_bins[i]-1;
			}
			cntMap[key]++;// rely on the default initialization
		}
		// convert count to frequency
		for (auto it=cntMap.begin(); it!=cntMap.end(); ++it) {
			pdf[it->first] = it->second*1.0/nElements;	
		}
		_computeJointCdf();

		// 2. update the mean 
		mean = ublas::zero_vector<Real>(num_comps);
		for (auto it=pdf.begin(); it!=pdf.end(); ++it) {
			std::vector<int> bin = it->first;
			for(int i=0; i<bin.size(); i++) {
				mean[i] += (min_vals[i]+(bin[i]+0.5)*bin_widths[i])*it->second; 
			}
		}

		// 3. update the cov matrix
		setCovariance( ublas::identity_matrix<Real>(num_comps) );
		for(int i=0; i<num_comps; i++) {
			for(int j=0; j<num_comps; j++) {
				for(int k=0; k<nElements; k++) {
					cov(i,j) += (dataAry[i][k]-mean[i])*(dataAry[j][k]-mean[j]);	
				}
			}
		}
		cov *= 1.0/nElements;
		setCovariance(cov);
	}

	__host__ __device__
	void _computeJointCdf() {
		assert(pdf.size()>0); 
		if(cdf.size()==pdf.size()) return; // already constructed before
		
		// sort PDF in decreasing order, i.e. move high frequency bins to front
		std::vector<std::pair<std::vector<int>, Real>> sorted_pdf(pdf.begin(), pdf.end());
		std::sort(sorted_pdf.begin(), sorted_pdf.end(), 
			boost::bind(&std::pair<std::vector<int>, Real>::second, _1) > 
			boost::bind(&std::pair<std::vector<int>, Real>::second, _2));
	
		// construct a CDF on the sparse JointHistogram
		Real accFreq = sorted_pdf[0].second;
		cdf.push_back(sorted_pdf[0]);	
		for(int i=1; i<sorted_pdf.size(); i++) {
			accFreq += sorted_pdf[i].second;
			std::pair<std::vector<int>, Real> tmp(sorted_pdf[i].first, accFreq);
			cdf.push_back(tmp);
		}
		const float epsilon = std::numeric_limits<float>::epsilon();	
		assert(abs(accFreq-1)<epsilon);
	}
};

// ------------------------------------------------------------------------------
// Below defines JointHistogram related generic functions
///
/// \brief Return joint mean of the joint histogram
///
__host__ __device__
inline std::vector<Real> getJointMean(const JointHistogram &dist) {
	std::vector<Real> m(dist.mean.size());
	std::copy(dist.mean.begin(), dist.mean.end(), m.begin());
	return m;
}

///
/// \brief Return PDF of x
///
///TODO: need to change the function name, getJointFrequency????
__host__ __device__
inline double getJointPdf(const JointHistogram &dist, const std::vector<Real> &x_) { 
	int comps = dist.getNumComps();
	std::vector<Real> mins = dist.getMinVals();
	std::vector<Real> widths = dist.getBinWidths();
	std::vector<int> bins = dist.getNumBins();
	std::vector<int> key(comps);
	for(int i=0; i<comps; i++) {// value to bin index
		key[i] = std::floor((x_[i]-mins[i])/widths[i]);
		if(key[i]>=0 && key[i]<bins[i])
			continue;	
		else
			key[i] = bins[i]-1;
	}
	return (dist.getDistr())[key]; 
}

///
/// \brief binary search vec and reture the index
///	
__host__
inline int biSearchNextLarge(std::vector<std::pair<std::vector<int>, Real>>& vec, Real p) {
	assert(vec.size()>0);
	if(vec.size()==1) return 0;
	// vec is sorted, increasing order
	int low = 0;
	int high = vec.size()-1;

	while(low<=high) {
		int mid = low+ (high-low)/2;
		if(vec[mid].second<p) low = mid+1;
		else if(vec[mid].second>p) high = mid-1;
		else return mid+1;
	}

	if(high<0) return 0; //p<vec[0];
	else{
		if( low > (vec.size()-1))
			return vec.size()-1; // p>=vec[vec.size()-1]
		else
			return (low<high) ? low+1 : high+1;
	}
}

///
/// \brief Return a random sample
///
__host__
inline std::vector<Real> getJointSample(const JointHistogram &dist) {
	if(dist.getDistr().size()==0) {
		printf("ERROR: distribution size is 0.\n");
		return std::vector<Real>(dist.getNumComps(), 0);
	}
	std::vector<std::pair<std::vector<int>, Real>> cdf;
	cdf = dist.getJointCdf();
	Real x = static_cast<Real>(rand())/RAND_MAX;
	int i=0;
	bool useBiSearch = true;
	if(useBiSearch) {
		i = biSearchNextLarge(cdf, x);
	}
	else {
		while(x > cdf[i].second)	
			i++;
	}

	// if reach here, we find the bin
	std::vector<int> bin = cdf[i].first;
	std::vector<Real> smp;
	std::vector<Real> mins = dist.getMinVals();
	std::vector<Real> widths = dist.getBinWidths();
	bool binCenter = true;
	for(int i=0; i<bin.size(); i++) { 
		Real smpPos;
		if(binCenter) smpPos = 0.5;// use bin center
		else smpPos = static_cast<Real>(rand())/RAND_MAX;
		smp.push_back(mins[i]+(bin[i]+smpPos)*widths[i]);
	}
	return smp;
}

///
/// \brief Return a random sample using random engine
///
__host__ __device__
inline std::vector<Real> getJointSample(const JointHistogram &dist, thrust::default_random_engine &rng) {
  // TODO version used on device side
	return getJointSample(dist);
}

///
/// \brief Print itself
///
__host__
inline std::ostream& operator<<(std::ostream& os, const JointHistogram &dist) {
    os <<  "<JointHistogram: mean=" << dist.mean << ", covariance=" << dist.getCovariance() << ">" ;
    return os;
}

__host__ __device__
inline std::string getName(const JointHistogram &dist) {
    return "JointHistogram";
}

///
/// \brief construct JointHistogram outside the class
///
__host__ __device__
inline JointHistogram eddaComputeJointHistogram(std::vector<Real*>& dataAry, int nElements, const std::vector<Real> &mins, 
				const std::vector<Real> &maxs, const std::vector<int> &nBins) {
	JointHistogram tmp(dataAry, nElements, mins, maxs, nBins); 
	return tmp;
}

}  // namespace dist
}  // namespace edda

#endif  // DIST_HISTOGRAM_H_
