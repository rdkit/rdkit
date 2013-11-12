/*******************************************************************************
shapeAlignment.cpp - Shape-it
 
Copyright 2012 by Silicos-it, a division of Imacosi BVBA
 
This file is part of Shape-it.

	Shape-it is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as published 
	by the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	Shape-it is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public License
	along with Shape-it.  If not, see <http://www.gnu.org/licenses/>.

Shape-it is linked against OpenBabel version 2.

	OpenBabel is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation version 2 of the License.

***********************************************************************/



#include <Shape/shapeAlignment.h>



ShapeAlignment::ShapeAlignment(GaussianVolume& gRef, GaussianVolume& gDb)
: _gRef(&gRef)
, _gDb(&gDb)
, _rAtoms(0)
, _rGauss(0)
, _dAtoms(0)
, _dGauss(0)
, _maxSize(0)
, _maxIter(50)
, _matrixMap()
{
	// Loop over the single atom volumes of both molecules and make the combinations
	_rAtoms = gRef.levels[0];
	_dAtoms = gDb.levels[0];
	_rGauss = gRef.gaussians.size();
	_dGauss = gDb.gaussians.size();
	_maxSize = _rGauss * _dGauss + 1;
}



ShapeAlignment::~ShapeAlignment(void)
{
	_gRef = NULL;
	_gDb = NULL;
	
	// Clear the matrix map
	for (MatIter mi = _matrixMap.begin(); mi != _matrixMap.end(); ++mi)
   {
		if (mi->second != NULL)
      {
			delete mi->second;
         mi->second = NULL;
		}
	}	
}



AlignmentInfo 
ShapeAlignment::gradientAscent(SiMath::Vector rotor)
{
	// Create a queue to hold the pairs to process
	std::queue<std::pair<unsigned int, unsigned int> > processQueue;
	
	// Helper variables
	// Overlap matrix is stored as a double error
	double* Aij;
	double Aq[4];
	
	double Vij(0.0);
   double qAq(0.0);
	
	std::vector<unsigned int> * d1(NULL);
	std::vector<unsigned int> * d2(NULL);
	std::vector<unsigned int>::iterator it1;
	
	// Gradient information
	SiMath::Vector overGrad(4, 0.0);
   
   // Hessian
	SiMath::Matrix overHessian(4, 4, 0.0);
	
	// Overlap volume
	double atomOverlap(0.0);
	double pharmOverlap(0.0);
		
	// Solution info
	AlignmentInfo res;
	
	double oldVolume(0.0);
	unsigned int iterations(0);
   unsigned int mapIndex(0);

//	unsigned int DBSize = _gDb->pharmacophores.size();

	while (iterations < 20)
   {
		// reset volume
		atomOverlap = 0.0;
		pharmOverlap = 0.0;
		iterations++;
				
		// reset gradient
		overGrad = 0.0;
      
		// reset hessian
		overHessian = 0.0;
				
		double lambda(0.0);
				
		// iterator over the matrix map
		MatIter matIter;
      
		// create atom-atom overlaps 
		for (unsigned int i(0); i < _rAtoms; ++i)
      {
			for (unsigned int j(0); j < _dAtoms; ++j)
         {
				mapIndex = (i * _dGauss) + j;
				
				if ((matIter = _matrixMap.find(mapIndex)) == _matrixMap.end())
            {
					Aij = _updateMatrixMap(_gRef->gaussians[i], _gDb->gaussians[j]);
               
					// add to map
					_matrixMap[mapIndex] = Aij;
				}
            else
            {
					Aij = matIter->second;
				}
				
				// rotor product
				Aq[0] =  Aij[0] * rotor[0] +  Aij[1] * rotor[1] +  Aij[2] * rotor[2] +  Aij[3] * rotor[3];
				Aq[1] =  Aij[4] * rotor[0] +  Aij[5] * rotor[1] +  Aij[6] * rotor[2] +  Aij[7] * rotor[3];
				Aq[2] =  Aij[8] * rotor[0] +  Aij[9] * rotor[1] + Aij[10] * rotor[2] + Aij[11] * rotor[3];
				Aq[3] = Aij[12] * rotor[0] + Aij[13] * rotor[1] + Aij[14] * rotor[2] + Aij[15] * rotor[3];
				
				qAq = rotor[0] * Aq[0] + rotor[1]*Aq[1] + rotor[2]*Aq[2] + rotor[3]*Aq[3];
				
				// compute overlap volume
				Vij = Aij[16] * exp( -qAq );
				
				// check if overlap is sufficient enough, should be more than 0.1 for atom - atom overlap
				if ( Vij/(_gRef->gaussians[i].volume + _gDb->gaussians[j].volume - Vij ) < EPS )
            {
					continue;
            }
				
				// add to overlap volume
				atomOverlap += Vij;

				// update gradient -2Vij (Aijq);
				double v2 = 2.0 * Vij;
				
				lambda -= v2 * qAq;
				
				overGrad[0] -= v2 * Aq[0];
				overGrad[1] -= v2 * Aq[1];
				overGrad[2] -= v2 * Aq[2];
				overGrad[3] -= v2 * Aq[3];
				
				// overHessian += 2*Vij(2*Aijq'qAij-Aij); (only upper triangular part)
				overHessian[0][0] += v2 * (2.0 * Aq[0]*Aq[0] - Aij[0]);
				overHessian[0][1] += v2 * (2.0 * Aq[0]*Aq[1] - Aij[1]);
				overHessian[0][2] += v2 * (2.0 * Aq[0]*Aq[2] - Aij[2]);
				overHessian[0][3] += v2 * (2.0 * Aq[0]*Aq[3] - Aij[3]);
				overHessian[1][1] += v2 * (2.0 * Aq[1]*Aq[1] - Aij[5]);
				overHessian[1][2] += v2 * (2.0 * Aq[1]*Aq[2] - Aij[6]);
				overHessian[1][3] += v2 * (2.0 * Aq[1]*Aq[3] - Aij[7]);
				overHessian[2][2] += v2 * (2.0 * Aq[2]*Aq[2] - Aij[10]);
				overHessian[2][3] += v2 * (2.0 * Aq[2]*Aq[3] - Aij[11]);
				overHessian[3][3] += v2 * (2.0 * Aq[3]*Aq[3] - Aij[15]);
				
				// loop over child nodes and add to queue
				d1 = _gRef->childOverlaps[i];
				d2 = _gDb->childOverlaps[j];
				
				// first add (i,child(j))
				if ( d2 != NULL )
            {
					for ( it1 = d2->begin(); it1 != d2->end(); ++it1 )
               {
						processQueue.push(std::make_pair<unsigned int, unsigned int>(i, *it1));
					}
				}
            
				// second add (child(i),j)
				if ( d1 != NULL )
            {
					for ( it1 = d1->begin(); it1 != d1->end(); ++it1 )
               {
						processQueue.push(std::make_pair<unsigned int, unsigned int>(*it1, j));
					}
				}
			}
      }

      while( !processQueue.empty() )
      {
         // get next element from queue
         std::pair<unsigned int, unsigned int> nextPair = processQueue.front();
         processQueue.pop();
					
         unsigned int i = nextPair.first;
         unsigned int j = nextPair.second;
					
         // check cache
         mapIndex = (i*_dGauss) + j;
         if ( (matIter = _matrixMap.find(mapIndex)) == _matrixMap.end() )
         {
            Aij = _updateMatrixMap(_gRef->gaussians[i], _gDb->gaussians[j]);
            _matrixMap[mapIndex] = Aij;
         }
         else
         {
            Aij = matIter->second;
         }
					
         // rotor product
         Aq[0] =  Aij[0] * rotor[0] +  Aij[1] * rotor[1] +  Aij[2] * rotor[2] +  Aij[3] * rotor[3];
         Aq[1] =  Aij[4] * rotor[0] +  Aij[5] * rotor[1] +  Aij[6] * rotor[2] +  Aij[7] * rotor[3];
         Aq[2] =  Aij[8] * rotor[0] +  Aij[9] * rotor[1] + Aij[10] * rotor[2] + Aij[11] * rotor[3];
         Aq[3] = Aij[12] * rotor[0] + Aij[13] * rotor[1] + Aij[14] * rotor[2] + Aij[15] * rotor[3];
					
         qAq = rotor[0] * Aq[0] + rotor[1]*Aq[1] + rotor[2]*Aq[2] + rotor[3]*Aq[3];
					
         // compute overlap volume
         Vij = Aij[16] * exp( -qAq );
					
         // check if overlap is sufficient enough
         if ( fabs(Vij)/(_gRef->gaussians[i].volume + _gDb->gaussians[j].volume-fabs(Vij)) < EPS )
         {
            continue;
         }
					
         atomOverlap += Vij;
					
         double v2 = 2.0*Vij;				
					
         lambda -= v2 * qAq;
					
         // update gradient -2Vij (Cij*Aij*q);
         overGrad[0] -= v2 * Aq[0];
         overGrad[1] -= v2 * Aq[1];
         overGrad[2] -= v2 * Aq[2];
         overGrad[3] -= v2 * Aq[3];
         
         // hessian 2*Vij(2*AijqqAij-Aij); (only upper triangular part)
         overHessian[0][0] += v2 * (2.0 * Aq[0]*Aq[0] - Aij[0]);
         overHessian[0][1] += v2 * (2.0 * Aq[0]*Aq[1] - Aij[1]);
         overHessian[0][2] += v2 * (2.0 * Aq[0]*Aq[2] - Aij[2]);
         overHessian[0][3] += v2 * (2.0 * Aq[0]*Aq[3] - Aij[3]);
         overHessian[1][1] += v2 * (2.0 * Aq[1]*Aq[1] - Aij[5]);
         overHessian[1][2] += v2 * (2.0 * Aq[1]*Aq[2] - Aij[6]);
         overHessian[1][3] += v2 * (2.0 * Aq[1]*Aq[3] - Aij[7]);
         overHessian[2][2] += v2 * (2.0 * Aq[2]*Aq[2] - Aij[10]);
         overHessian[2][3] += v2 * (2.0 * Aq[2]*Aq[3] - Aij[11]);
         overHessian[3][3] += v2 * (2.0 * Aq[3]*Aq[3] - Aij[15]);
					
					
         // loop over child nodes and add to queue
         d1 = _gRef->childOverlaps[i];
         d2 = _gDb->childOverlaps[j];
         if ( d1 != NULL && _gRef->gaussians[i].nbr >  _gDb->gaussians[j].nbr )
         {
            for ( it1 = d1->begin(); it1 != d1->end(); ++it1 )
            {
               // add (child(i),j)
               processQueue.push(std::make_pair<unsigned int, unsigned int>(*it1, j));
            }
         }
         else
         {
            // first add (i,child(j))
            if ( d2 != NULL )
            {
               for ( it1 = d2->begin(); it1 != d2->end(); ++it1 )
               {
                  processQueue.push(std::make_pair<unsigned int, unsigned int>(i, *it1));
               }
            }
            if ( d1 != NULL && _gDb->gaussians[j].nbr - _gRef->gaussians[i].nbr < 2 )
            {
               for ( it1 = d1->begin(); it1 != d1->end(); ++it1 )
               {
                  // add (child(i),j)
                  processQueue.push(std::make_pair<unsigned int, unsigned int>(*it1, j));
               }
            }
         }
      }
				

      // check if the new volume is better than the previously found one
      // if not quit the loop
      if ( iterations > 6 && atomOverlap < oldVolume + 0.0001 )
      {
         break;
      }
				
      // store latest volume found
      oldVolume = atomOverlap;
				
      // no measurable overlap between two volumes
      if ( std::isnan(lambda) || std::isnan(oldVolume) || oldVolume == 0 )
      {
         break;
      }
				
      // update solution 
      if ( oldVolume > res.overlap)
      {
         res.overlap = atomOverlap;
         res.rotor = rotor;
         if ( res.overlap/(_gRef->overlap + _gDb->overlap - res.overlap ) > 0.99 )
         {
            break;
         }
      }
				
      // update the gradient and hessian 
      overHessian -= lambda;
				
      // fill lower triangular of the hessian matrix
      overHessian[1][0] = overHessian[0][1];
      overHessian[2][0] = overHessian[0][2];
      overHessian[2][1] = overHessian[1][2];
      overHessian[3][0] = overHessian[0][3];
      overHessian[3][1] = overHessian[1][3];
      overHessian[3][2] = overHessian[2][3];
				
      // update gradient to make h
      overGrad[0] -= lambda * rotor[0];
      overGrad[1] -= lambda * rotor[1];
      overGrad[2] -= lambda * rotor[2];
      overGrad[3] -= lambda * rotor[3];
      
      // update gradient based on inverse hessian
      SiMath::Vector tmp = rowProduct(overHessian, overGrad);
      double h = overGrad[0]*tmp[0] + overGrad[1]*tmp[1] + overGrad[2]*tmp[2] + overGrad[3]*tmp[3];
      
      // small scaling of the gradient
      h = 1/h*(overGrad[0]*overGrad[0] + overGrad[1]*overGrad[1] + overGrad[2]*overGrad[2] + overGrad[3]*overGrad[3]);
      overGrad *= h;
				
      // update rotor based on gradient information
      rotor -= overGrad;
				
      // normalise rotor such that it has unit norm
      double nr = sqrt(rotor[0]*rotor[0] + rotor[1]*rotor[1] + rotor[2]*rotor[2] + rotor[3]*rotor[3]);
				
      rotor[0] /= nr;
      rotor[1] /= nr;
      rotor[2] /= nr;
      rotor[3] /= nr;
      
	} // end of endless while loop
		return res;
}



AlignmentInfo 
ShapeAlignment::simulatedAnnealing(SiMath::Vector rotor)
{
	// create a queue to hold the pairs to process
	std::queue<std::pair<unsigned int, unsigned int> > processQueue;
	
	// helper variables
	// overlap matrix is stored as a double error
	double* Aij;
	double Aq[4];
	
	// map store store already computed matrices
	MatIter matIter;
	
	double Vij(0.0), qAq(0.0);
	
	std::vector<unsigned int>* d1(NULL);
	std::vector<unsigned int>* d2(NULL);
	std::vector<unsigned int>::iterator it1;
	
	// overlap volume
	double atomOverlap(0.0);
	double pharmOverlap(0.0);
	
	double dTemperature = 1.1;
	
	// solution info
	AlignmentInfo res;
	
	SiMath::Vector oldRotor(rotor);
	SiMath::Vector bestRotor(rotor);
		
	double oldVolume(0.0);
	double bestVolume(0.0);
	unsigned int iterations(0);
	unsigned int sameCount(0);
	unsigned int mapIndex(0);
//	unsigned int DBSize = _gDb->pharmacophores.size();
		
	while ( iterations < _maxIter)
   {
		// reset volume
		atomOverlap = 0.0;
		pharmOverlap = 0.0;
		
		++iterations;

		// temperature of the simulated annealing step
		double T = sqrt((1.0+iterations)/dTemperature);
		
		// create atom-atom overlaps 
		for ( unsigned int i=0; i<_rAtoms; ++i )
      {
			for ( unsigned int j=0; j<_dAtoms; ++j )
         {
				mapIndex = (i * _dGauss) + j;
				
				if ( (matIter = _matrixMap.find(mapIndex)) == _matrixMap.end() )
            {
					Aij = _updateMatrixMap(_gRef->gaussians[i], _gDb->gaussians[j]);
					_matrixMap[mapIndex] = Aij;
				}
            else
            {
					Aij = matIter->second;
				}
				
				// rotor product
				Aq[0] =  Aij[0] * rotor[0] +  Aij[1] * rotor[1] +  Aij[2] * rotor[2] +  Aij[3] * rotor[3];
				Aq[1] =  Aij[4] * rotor[0] +  Aij[5] * rotor[1] +  Aij[6] * rotor[2] +  Aij[7] * rotor[3];
				Aq[2] =  Aij[8] * rotor[0] +  Aij[9] * rotor[1] + Aij[10] * rotor[2] + Aij[11] * rotor[3];
				Aq[3] = Aij[12] * rotor[0] + Aij[13] * rotor[1] + Aij[14] * rotor[2] + Aij[15] * rotor[3];
				
				qAq = rotor[0] * Aq[0] + rotor[1]*Aq[1] + rotor[2]*Aq[2] + rotor[3]*Aq[3];
				
				// compute overlap volume
				Vij = Aij[16] * exp( -qAq );
				
				// check if overlap is sufficient enough, should be more than 0.1 for atom - atom overlap
				if ( Vij/(_gRef->gaussians[i].volume + _gDb->gaussians[j].volume - Vij ) < EPS )
            {
					continue;
            }
				
				// add to overlap volume
				atomOverlap += Vij;
				
				// loop over child nodes and add to queue
				d1 = _gRef->childOverlaps[i];
				d2 = _gDb->childOverlaps[j];
				
				// first add (i,child(j))
				if ( d2 != NULL )
            {
					for ( it1 = d2->begin(); it1 != d2->end(); ++it1 )
               {
						processQueue.push(std::make_pair<unsigned int, unsigned int>(i, *it1));
					}
				}
				// second add (child(i),j)
				if ( d1 != NULL )
            {
					for ( it1 = d1->begin(); it1 != d1->end(); ++it1 )
               {
						processQueue.push(std::make_pair<unsigned int, unsigned int>(*it1, j));
					}
				}
			}
		}

		while( !processQueue.empty() )
		{
			// get next element from queue
			std::pair<unsigned int, unsigned int> nextPair = processQueue.front();
			processQueue.pop();
			
			unsigned int i = nextPair.first;
			unsigned int j = nextPair.second;
			
			// check cache
			mapIndex = (i*_dGauss) + j;
			if ( (matIter = _matrixMap.find(mapIndex)) == _matrixMap.end() )
         {
				Aij = _updateMatrixMap(_gRef->gaussians[i],_gDb->gaussians[j]);
				_matrixMap[mapIndex] = Aij;
			}
         else
         {
				Aij = matIter->second;
			}
			
			// rotor product
			Aq[0] =  Aij[0] * rotor[0] +  Aij[1] * rotor[1] +  Aij[2] * rotor[2] +  Aij[3] * rotor[3];
			Aq[1] =  Aij[4] * rotor[0] +  Aij[5] * rotor[1] +  Aij[6] * rotor[2] +  Aij[7] * rotor[3];
			Aq[2] =  Aij[8] * rotor[0] +  Aij[9] * rotor[1] + Aij[10] * rotor[2] + Aij[11] * rotor[3];
			Aq[3] = Aij[12] * rotor[0] + Aij[13] * rotor[1] + Aij[14] * rotor[2] + Aij[15] * rotor[3];
			
			qAq = rotor[0] * Aq[0] + rotor[1]*Aq[1] + rotor[2]*Aq[2] + rotor[3]*Aq[3];
			
			// compute overlap volume
			Vij = Aij[16] * exp( -qAq );
			
			// check if overlap is sufficient enough
			if ( fabs(Vij)/(_gRef->gaussians[i].volume + _gDb->gaussians[j].volume-fabs(Vij)) < EPS ) 
         {
				continue;
         }
			
			// even number of overlap atoms => addition to volume
			atomOverlap += Vij;
			
			// loop over child nodes and add to queue
			d1 = _gRef->childOverlaps[i];
			d2 = _gDb->childOverlaps[j];
			if ( d1 != NULL && _gRef->gaussians[i].nbr >  _gDb->gaussians[j].nbr )
         {
				for ( it1 = d1->begin(); it1 != d1->end(); ++it1 )
            {
					// add (child(i),j)
					processQueue.push(std::make_pair<unsigned int, unsigned int>(*it1, j));
				}
			}
         else
         {
				// first add (i,child(j))
				if ( d2 != NULL )
            {
					for ( it1 = d2->begin(); it1 != d2->end(); ++it1 )
               {
						processQueue.push(std::make_pair<unsigned int, unsigned int>(i, *it1));
					}
				}
				if ( d1 != NULL && _gDb->gaussians[j].nbr - _gRef->gaussians[i].nbr < 2 ) 
            {
					for ( it1 = d1->begin(); it1 != d1->end(); ++it1 )
               {
						// add (child(i),j)
						processQueue.push(std::make_pair<unsigned int, unsigned int>(*it1, j));
					}
				}
			}
		}

		// check if the new volume is better than the previously found one
		double overlapVol = atomOverlap;
		if ( overlapVol < oldVolume )
		{
			double D = exp(-sqrt(oldVolume - overlapVol))/T;
			if ( SiMath::randD(0,1) < D )
         {
				oldRotor = rotor;
				oldVolume = overlapVol;
				sameCount = 0;
			}
         else
         {
				++sameCount;
				if ( sameCount == 30 )
            {
					iterations = _maxIter;
            }
			}
		}
      else
      {
			// store latest volume found
			oldVolume = overlapVol;
			oldRotor = rotor;
         
			// update best found so far
			bestRotor = rotor;
			bestVolume = overlapVol;
			sameCount = 0;
			
			// check if it is better than the best solution found so far
			if (overlapVol > res.overlap)
         {
				res.overlap = atomOverlap;
				res.rotor = rotor;
				if ((res.overlap / (_gRef->overlap + _gDb->overlap - res.overlap )) > 0.99)
            {
					break;
            }
			}
			
		}
		
		// make random permutation
		// double range = 0.05;
		double range = 0.1/T;
		rotor[0] = oldRotor[0] + SiMath::randD(-range, range);
		rotor[1] = oldRotor[1] + SiMath::randD(-range, range);
		rotor[2] = oldRotor[2] + SiMath::randD(-range, range);
		rotor[3] = oldRotor[3] + SiMath::randD(-range, range);
		
		// normalise rotor such that it has unit norm
		double nr = sqrt(rotor[0]*rotor[0] + rotor[1]*rotor[1] + rotor[2]*rotor[2] + rotor[3]*rotor[3]);
		rotor[0] /= nr;
		rotor[1] /= nr;
		rotor[2] /= nr;
		rotor[3] /= nr;
		
	} // end of endless while loop

	return res;
}



void
ShapeAlignment::setMaxIterations(unsigned int i)
{
   _maxIter = i;
   return;
}



double* 
ShapeAlignment::_updateMatrixMap(AtomGaussian& a, AtomGaussian& b)
{
	double * A = new double[17];
	
	// variables to store sum and difference of components
	double dx = (a.center.x - b.center.x);
	double dx2 = dx * dx;
	double dy = (a.center.y - b.center.y);
	double dy2 = dy * dy;
	double dz = (a.center.z - b.center.z);
	double dz2 = dz * dz;
	double sx = (a.center.x + b.center.x);
	double sx2 = sx * sx;
	double sy = (a.center.y + b.center.y);
	double sy2 = sy * sy;
	double sz = (a.center.z + b.center.z);
	double sz2 = sz * sz;
	
	// update overlap matrix
	double C = a.alpha * b.alpha /(a.alpha + b.alpha);
	A[0] = C * (dx2 + dy2 + dz2);
	A[1] = C * (dy*sz - dz*sy);
	A[2] = C * (dz*sx - dx*sz);
	A[3] = C * (dx*sy - dy*sx);
	A[4] = A[1];
	A[5] = C * (dx2 + sy2 + sz2);
	A[6] = C * (dx*dy - sx*sy);
	A[7] = C * (dx*dz - sx*sz);
	A[8] = A[2];
	A[9] = A[6];
	A[10] = C * (sx2 + dy2 + sz2);
	A[11] = C * (dy*dz - sy*sz);
	A[12] = A[3];
	A[13] = A[7];
	A[14] = A[11];
	A[15] = C * (sx2 + sy2 + dz2);
   
	// last elements holds overlap scaling constant
	// even number of overlap atoms => addition to volume
	// odd number => substraction
	if ( ( a.nbr + b.nbr ) % 2 == 0 )
   {
		A[16] = a.C * b.C * pow(PI/(a.alpha + b.alpha),1.5);
	}
   else
   {
		A[16] = -a.C * b.C * pow(PI/(a.alpha + b.alpha),1.5);						
	}
	
	return A;
}
