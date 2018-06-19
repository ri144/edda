// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef DIST_JOINT_GAUSSIAN_H_
#define DIST_JOINT_GAUSSIAN_H_

#include <cstdlib>
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES  // For Visual Studio
#include <math.h>

#include <boost/limits.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>
#include <boost/random.hpp>

#include "common.h"
#include "distribution_tag.h"
#include "core/statistics.h"
#include "invert_matrix.h"

#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

namespace edda {
	namespace dist {

		// ------------------------------------------------------------------------------
		///
		/// \brief Defines a Gaussian class
		///
		struct EDDA_EXPORT JointGaussian : public ContinuousDistributionTag, public JointDistributionTag {

			// constructor
			__host__ __device__
				JointGaussian() {
					mean = ublas::zero_vector<Real>(3);
					setMatrices(ublas::identity_matrix<Real>(3));
				}

			///
			/// \brief Constructor
			///
			__host__ __device__
				JointGaussian(const ublas_vector &mean, const ublas_matrix &cov) {
					assert(mean.size() == cov.size1() && mean.size() == cov.size2());
					this->mean = mean;
					setMatrices(cov);
				}


                        ///
                        /// \brief computer lower matrix and precision cholesky matrix for probability estimat and sampling
                        /// \param cov the input covariance matrix
                        ///
                        __host__ __device__
                                void setMatrices(const ublas_matrix &cov) {
                                        //Decomposition computation steps
                                        //1. use cholescky decomposition: covMat -> L (lower matrix)
                                        //2. solve LX = I to get X (this can be solve by Lx=e column by column)
                                        //3. put transpose(X) in precCholMat
                                        //4. calculate logPDet from precCholMat
                                        //Note: step 1 and 2 now is on Eigen libray
                                        //      they are possible to be implement by ourself, 
                                        //      then we can get rid of Eigen library.
                                        this->cov = cov;  //keep here

                                        int covSize = this->cov.size1(); //cov.size1 == cov.size2

                                        //convert matrix data structure
                                        MatrixXd rawCov(covSize, covSize);
                                        for (int i = 0; i < covSize; i++){
                                                for (int j = 0; j < covSize; j++){
                                                        rawCov(i,j) = this->cov(i, j);
                                                }
                                        }

                                        //Chol decomposition and triangular solver from Eigen
                                        //for probability computation (EM Need this)
                                        Eigen::LLT<Eigen::MatrixXd>  chol(rawCov);
                                        //Handle non positive definite covariance somehow:
                                        if(chol.info()!=Eigen::Success){
                                            std::cout << "Choleschy Decomposition Error!! (use double precision may solve this problem)" << std::endl;
                                            exit(0);
                                        }
                                        MatrixXd L = chol.matrixL();
                                        lowerMat.resize(covSize, covSize);
                                        for (int i = 0; i < covSize; i++){
                                                for (int j = 0; j < covSize; j++){
                                                        lowerMat(i,j) = L(i, j);
                                                }
                                        }
                                        MatrixXd I = MatrixXd::Identity(covSize,covSize);
                                        MatrixXd preciChol = L.triangularView<Lower>().solve(I);
                                        precCholMat.resize(covSize, covSize); 
                                        for (int i = 0; i < covSize; i++){
                                            for (int j = 0; j < covSize; j++){
                                                precCholMat(i,j) = preciChol(j,i);     
                                            }
                                        }

                                        //Compute logPDet
                                        logPDet = 0;
                                        for( int i=0; i<covSize; i++ ){
                                            logPDet += log( precCholMat(i,i) ); //diag elements
                                        }
                                }

            
			///
			/// \brief Return a sample drawn from this joint Gaussian
			/// \param rng random engine
			///
			__host__ __device__
				std::vector<Real> getJointSample(/*thrust::default_random_engine &rng*/) const {
					//draw nVar samples from independent standard normal variant
					//thrust::random::normal_distribution<Real> ndist(0,1);
					ublas_vector r(cov.size1());

					for (int j = 0; j < r.size(); j++){
						r(j) = static_cast<Real>(rand())/RAND_MAX;
					}

					//transform the above sample by the covariance matrix
					ublas_vector s(cov.size1());
					boost::numeric::ublas::axpy_prod(lowerMat, r, s, true);

					std::vector<Real> retS(cov.size1());
					for (int j = 0; j < r.size(); j++){
						retS[j] = s(j) + mean(j);
					}

					return retS;
				}

			///
			/// \brief Return log probability of x (This function is used by EM)
			/// \param x_ vector of a sample for probability estimation
			///
			__host__ __device__
				inline Real getJointLogPdf(const std::vector<Real> x_) const
			{
					int k = this->mean.size();
					assert(x_.size() == k);
					ublas_vector x(k);
					std::copy(x_.begin(), x_.end(), x.begin());
					x = x - this->mean;					
					ublas_vector tmp = ublas::prod(ublas::trans(x), getPrecCholMat());
					Real density = ublas::inner_prod(tmp, tmp);

					return -0.5 * (k * 1.83787706641 + density ) + getLogDet();
			}

			__host__ __device__
				inline Real getJointPdf(const std::vector<Real> x_)
			{
				int k = this->mean.size();
				assert(x_.size() == k);
				ublas_vector x(k);
				std::copy(x_.begin(), x_.end(), x.begin());
				x = x - this->mean;
				ublas_vector tmp = ublas::prod(ublas::trans(x), getPrecCholMat());
				Real density = ublas::inner_prod(tmp, tmp);

				return exp( -0.5 * (k * 1.83787706641 + density ) + getLogDet());
			}

			///
			/// \brief Return mean vector of this Gaussian
			///
			__host__ __device__
				const ublas_vector &getMean() const { return this->mean; }

			///
			/// \brief Return covariance matrix of this Gaussian
			///
			__host__ __device__
				const ublas_matrix &getCovariance() const { return this->cov; }

			///
			/// \brief Return eigen matrix of this Gaussian's covariance matrix
			///
			__host__ __device__
				const ublas_matrix &getPrecCholMat() const { return this->precCholMat; }

			///
			/// \brief Return U matrix of this Gaussian's covariance matrix
			///
			__host__ __device__
				Real getLogDet() const { return this->logPDet; }

		private:
			ublas_vector mean; 	    // boost's vector for mean vector
			ublas_matrix cov;	    //covariance matrix
			ublas_matrix precCholMat;   //precisionCholeskyMat:cov->L and solve LX=I, the X is this
                        ublas_matrix lowerMat;      //lower matrix from covMat cholescky decomposition
			Real logPDet;		    //log determinate of covariance matrix
		};

		// ------------------------------------------------------------------------------
		// Below defines JointGaussian related generic functions
		__host__ __device__
			inline std::vector<Real> getJointMean(const JointGaussian &dist)
		{
			ublas_vector mean = dist.getMean();
			std::vector<Real> m(mean.size());
			std::copy(mean.begin(), mean.end(), m.begin());
			return m;
		}

		///
		/// \brief Return PDF of x
		/// \param dist a distribution (Gaussian)
		/// \param x_ vector of a sample for probability estimation
		///
		__host__ __device__
			inline Real getJointPdf(const JointGaussian &dist, const std::vector<Real> x_)
		{
			ublas_vector mean = dist.getMean();
			int k = mean.size();
			assert(x_.size() == k);
			ublas_vector x(k);  
			std::copy(x_.begin(), x_.end(), x.begin());
			x = x - mean;
			ublas_vector tmp= ublas::prod(ublas::trans(x), dist.getPrecCholMat());
			Real density = ublas::inner_prod(tmp, tmp);

			return exp(-0.5 * (k * 1.83787706641 + density ) + dist.getLogDet());
		}

		///
		/// \brief Return a random sample using random engine
		/// \param dist a distribution (Gaussian)
		/// \param rng random engine
		///
		__host__ __device__
			inline std::vector<Real> getJointSample(const JointGaussian &dist/*, thrust::default_random_engine &rng*/)
		{
			//return getJointSample(dist, rng);
			return dist.getJointSample();
		}

		///
		/// \brief Print itself
		/// \param os outstream
		/// \param dist a distribution
		///
		__host__
			inline std::ostream& operator<<(std::ostream& os, const JointGaussian &dist)
		{
				os << "<JointGaussian: mean=" << dist.getMean() << ", covariance=" << dist.getCovariance() << ", >";
				return os;
		}

		///
		/// \brief Print the distribution type
		/// \param dist a distribution
		///
		__host__ __device__
		inline std::string getName(const JointGaussian &dist) {
			return "JointGaussian";
		}


	}  // namespace dist

	///
	/// \brief Create a Joint Gaussian
	/// \param dataAry input sample vectors for modeling
	/// \param nSamples number of input samples (vectors)
	///
	inline dist::JointGaussian eddaComputeJointGaussian(std::vector<Real*>& dataAry, int nSamples)
	{
		int nVar = dataAry.size();
		if( nSamples <= 1 ){
			std::cout<<"Number of samples must be larger than 1, return a non-initialized gaussian"<<std::endl;
        	return dist::JointGaussian();
		}

		ublas_vector mean(nVar);
		for (int j = 0; j < nVar; j++){
			mean(j) = 0;
			for (int i = 0; i < nSamples; i++)
				mean(j) += dataAry[j][i];
			mean(j) /= static_cast<Real>(nSamples);
		}

		ublas_matrix cov(nVar, nVar);
		for( int i=0; i<nVar; i++ ){
			for( int j=0; j<nVar; j++ ){
				cov(i,j) = 0;
				for( int k=0; k<nSamples; k++ ){
					cov(i,j) += (dataAry[i][k] - mean(i)) * (dataAry[j][k] - mean(j));
				}
				cov(i,j) /= static_cast<Real>(nSamples-1);
			}
		}

		return dist::JointGaussian(mean, cov);
	}
}  // namespace edda

#endif  // DIST_GAUSSIAN_H_
