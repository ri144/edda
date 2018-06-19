// Copyright 2015 The Edda Authors. All rights reserved.
// Use of this source code is governed by a MIT-style license that can be
// found in the LICENSE file.

#ifndef JointGMM_H
#define JointGMM_H

#include <iostream>
#include <typeinfo>
#include <iostream>

#include "common.h"
#include "edda_export.h"
#include "core/vector_matrix.h"
#include "core/tuple.h"
#include "core/thrust_common.h"
#include "distribution_tag.h"
#include "core/shared_ary.h"
#include "distributions/joint_gaussian.h"

#include "common.h"
#include "distribution_tag.h"
#include "core/statistics.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/assignment.hpp> 
//using namespace boost::numeric::ublas;

#define EMEPS 0.0000000000222044604

namespace edda {
	namespace dist {

		///
		/// \brief Define a Joint GMM class
		///
		class EDDA_EXPORT JointGMM : public ContinuousDistributionTag, public JointDistributionTag {		
		public:

			///
			/// \brief Constructor
			///
			JointGMM(){}

			///
			/// \brief Constructor by giving GMM's parameters
			/// \param gmm a given GMM
			///
			JointGMM(const JointGMM &gmm)
				: weights(gmm.weights),	gaus(gmm.gaus),	nVar(gmm.nVar),	nComp(gmm.nComp)
			{
			}

			///
			/// \brief Constructor by giving GMM's parameters and Gaussians
			/// \param _weights vector of wegiths
			/// \param _gaus gaussian vectors which contains multiple Gaussians
			/// \param _nVar number of variante
			/// \param _nComp number of Gaussian components
			///
			JointGMM(std::vector<Real> _weights, std::vector<JointGaussian> _gaus, int _nVar, int _nComp)
				: weights(_weights), gaus(_gaus), nVar(_nVar), nComp(_nComp)
			{
			}

			///
			/// \brief Set a Joint GMM
			/// \param _nVar number of variante
			/// \param _nComp number of Gaussian components
			/// \param _weights vector of wegiths
			/// \param _means mean vectors
			/// \param _covs covariance matrixs
			///
			void setGMM(int _nVar, int _nComp, ublas_vector &_weights, ublas_matrix &_means, ublas_matrix &_covs ){
				nVar = _nVar;
				nComp = _nComp;

				gaus.clear();
				weights.resize(nComp);
				std::copy(_weights.begin(), _weights.end(), weights.begin());

				for (int i = 0; i < nComp; i++){
					boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<Real> > m(_means, i);
					ublas_matrix cov_ = subrange( _covs, i*nVar, (i+1)*nVar, 0, nVar );

					JointGaussian g(m, cov_);
					gaus.push_back(g);
				}
			}

			///
			/// \brief Return means of a Gaussian component
			/// \param comp index of the Gaussian component
			///
			std::vector<Real> getMean(int comp) const{
				return getJointMean(gaus[comp]);
			}

			///
			/// \brief Return probability of x (sample)
			/// \param x samples (a vector)
			///
			Real getPdf(const std::vector<Real> x) const{
				Real pdf = 0;
				for (int i = 0; i < nComp; i++){
					pdf += weights[i] * getJointPdf(gaus[i], x);
				}

				return pdf;
			}

			///
			/// \brief Return weight probability of each component (only used by EM)
			/// \param x samples (a vector)
			///
			std::vector<Real> getCompWgtLogProbability(const std::vector<Real> x) const{
				std::vector<Real> wProbs(nComp, 0);
				for (int i = 0; i < nComp; i++){
					wProbs[i] = log(weights[i]) + gaus[i].getJointLogPdf(x); 
				}

				return wProbs;
			}

			///
			/// \brief Return a samplel drawn from GMM
			/// \param rng random engine
			///
			std::vector<Real> getJointSample(/*thrust::default_random_engine &rng*/) const{
				Real r = rand() / (Real)RAND_MAX;
				Real sum = 0;
				std::vector<Real> s;

				for (int i = 0; i < gaus.size(); i++){
					sum += weights[i];
					if (sum >= r || i==gaus.size()-1){
						s = gaus[i].getJointSample(/*rng*/);
						break;
					}
				}
				return s;
			}

			///
			/// \brief Print GMM parameters
			/// \param os outstream
			///
			void output(std::ostream& os) const{
				for (int i = 0; i < nComp; i++){
					os << "Component: " << i << std::endl;
					os << "weight: " << weights[i] << std::endl;
					os << gaus[i] << std::endl;
				}
			}

			///
			/// \brief Return number of variante of a Gaussian
			///
			int getNumVariables(){ return nVar; };

			///
			/// \brief Return number of components of a Gaussian
			///
			int getNumComponents(){ return nComp; };

			///
			/// \brief Return weight of components of a Gaussian
			///
			Real getWeight(int i){ return weights[i]; };

			///
			/// \brief Return one Gaussian of the GMM
			/// \param i the index of the return Gaussian component
			///
			JointGaussian getJointGaussian(int i) { return gaus[i]; };
		protected:
			std::vector<Real> weights;			// weight vector of GMM
			std::vector<JointGaussian> gaus;	// Gaussians of GMM
			int nVar;							// Number of variante
			int nComp;							// Number of Gaussian components
		};

		///
		/// \brief Return probability of x
		/// \param dist a distribution
		/// \param x a input sample
		///
		inline Real getPdf(const JointGMM &dist, const std::vector<Real> x)
		{
			return (Real)dist.getPdf(x);
		}

		///
		/// \brief Return a sample drawn from dist
		/// \param dist a distribution
		/// \param rng random engine
		///
		inline std::vector<Real> getJointSample(const JointGMM &dist/*, thrust::default_random_engine &rng*/)
		{
			return dist.getJointSample(/*rng*/);
		}

		///
		/// \brief Print itself
		/// \param os outstream
		/// \param dist a distribution
		///
		inline std::ostream& operator<<(std::ostream& os, const JointGMM &dist)
		{
			dist.output(os);
			return os;
		}

		///
		/// \brief Print the distribution typw
		/// \param x a distribution
		///
		__host__ __device__
		inline std::string getName(const JointGMM &x) {
			return "JointGMM";
		}

	}  // namespace dist

        void estimateParameters(ublas_matrix samples ,ublas_matrix resposibilities, int nSamples, int nComp, int nVar, 
                                ublas_vector& w, ublas_matrix& m, ublas_matrix& covs)
        {
                        for (int j = 0; j < nComp; j++){
                            w(j) = 0;
                            for (int smpPtr = 0; smpPtr < nSamples; smpPtr++){
                                    w(j) += resposibilities(smpPtr, j);
                            }
                            w(j) += 10 * std::numeric_limits<Real>::epsilon();
                        }


                        ublas_matrix weighted_X_sum = boost::numeric::ublas::prod(boost::numeric::ublas::trans(resposibilities) , samples);
                        for (int c = 0; c < nComp; c++){
                                subrange(m, c, c + 1, 0, nVar) = subrange(weighted_X_sum, c, c + 1, 0, nVar)  / w(c);
                        }

                        for (int c = 0; c < nComp; c++){
                                ublas_vector post(nSamples);
                                for (int s = 0; s < nSamples; s++){
                                        post(s) = resposibilities(s, c);
                                }
                                ublas_vector mu = row(m, c);
                                ublas_matrix diff(nSamples, nVar);
                                for (int s = 0; s < nSamples; s++){
                                        row(diff, s) = row(samples, s) - mu;
                                }

                                ublas_matrix postDiff( nVar, nSamples );
                                for (int s = 0; s < nSamples; s++){
                                        column(postDiff, s) = row(diff, s) * post(s);
                                }
                                ublas_matrix avg_cv = boost::numeric::ublas::prod(postDiff, diff);
                                avg_cv = avg_cv / w(c);
                                ublas_matrix cv = avg_cv + 0.000001 * ublas::identity_matrix<Real>(nVar);
                                subrange(covs, c*nVar, (c + 1)*nVar, 0, nVar) = cv;
                        }

                        for (int j = 0; j < nComp; j++){
                            w(j) /= (Real)nSamples;
                        }
        }//End estimateParameter()
	
	///
	/// \brief Compute the Gaussian Mixture Model 
	/// \param dataAry input sample vectors for modeling
	/// \param nSamples number of input samples (vectors)
	/// \param nComp number of Gaussian component of GMM
	///
	inline dist::JointGMM eddaComputeJointGMM(std::vector<Real*>& dataAry, int nSamples, int nComp, int verbose = 0, int n_iter=100, Real converge_threshold = 0.001)
	{
		int nVar = dataAry.size();

		//convert input sample to ublas_matrix
		ublas_matrix samples(nSamples, nVar);
		for (int j = 0; j < nVar; j++){
			for (int i = 0; i < nSamples; i++){
				samples(i, j) = dataAry[j][i];
			}
		}

		dist::JointGMM gmm;
		//initialization: 
		ublas_vector w(nComp);
		ublas_matrix m(nComp, nVar);
		ublas_matrix covs(nComp* nVar, nVar);
                ublas_matrix resposibilities(nSamples, nComp);
                Real current_total_log_norm = NAN;
                bool converged_ = false;

                //random resposibility to init weights means and covs
                for( int i=0; i<nSamples; i++ ){
                    Real sum = 0;
                    for( int j=0; j<nComp; j++ ){
                        int r = rand();
                        resposibilities(i,j) = r;
                        sum += resposibilities(i,j);
                    }
                    for(int j=0; j<nComp; j++ ) resposibilities(i,j) /= sum;
                }
                estimateParameters(samples, resposibilities, nSamples, nComp, nVar, w, m ,covs);

		for (int i = 0; i < n_iter; i++){
			gmm.setGMM(nVar, nComp, w, m, covs);
			Real prev_total_log_norm = current_total_log_norm;
                        if( verbose > 1 ){
                            std::cout << "Weights: " << w << std::endl;
                            std::cout << "Means: " << m << std::endl;
                            std::cout << "Covs: " << covs << std::endl << std::endl << std::endl;
                        }

			//E step (calculate resposibility)
			current_total_log_norm = 0;
			for (int smpPtr = 0; smpPtr < nSamples; smpPtr++){
				std::vector<Real> s(nVar);
				boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<Real> > mr(samples, smpPtr);
				std::copy(mr.begin(), mr.end(), s.begin());
				std::vector<Real> score = gmm.getCompWgtLogProbability(s);

                                Real log_prob_norm = 0;
                                for( int j = 0; j<nComp; j++ ){
                                    if( exp(score[j]) < std::numeric_limits<float>::min() ) {
                                        score[j] = log(std::numeric_limits<float>::min());
                                    }
                                    log_prob_norm += exp( score[j] );
                                }
                                log_prob_norm = log( log_prob_norm );
                                for( int j=0; j<nComp; j++ ){
                                    resposibilities( smpPtr, j ) = exp( score[j] - log_prob_norm );
                                }
                                current_total_log_norm += log_prob_norm;
			}

                        // Check converge
			current_total_log_norm /= (Real)nSamples;
			if (prev_total_log_norm == prev_total_log_norm){
				Real change = abs(current_total_log_norm - prev_total_log_norm);
				if (change < converge_threshold){//threshold for stopping EM (should be a user input para)
                                        converged_ = true;
					break;
				}
			}

			//M step (use resposibility to update weights, means, covs
                        estimateParameters(samples, resposibilities, nSamples, nComp, nVar, w, m ,covs);

                        if( verbose > 0 ){
                            std::cout << "Iters: " << i << ", log_norm: " << current_total_log_norm  << std::endl;
                        }
                        if( verbose > 1 ){
                            std::cout << "Weights: " << w << std::endl;
                            std::cout << "Means: " << m << std::endl;
                            std::cout << "Covs: " << covs << std::endl << std::endl << std::endl;
                        }
		}

		gmm.setGMM(nVar, nComp, w, m, covs);
		return gmm;
	}

}  // namespace edda

#endif // JointGMM_H
