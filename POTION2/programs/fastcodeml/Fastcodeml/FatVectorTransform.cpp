
#ifdef NEW_LIKELIHOOD

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cstdio>
#include <set>
#include "FatVectorTransform.h"
#include "MatrixSize.h"
#include "Exceptions.h"
#ifdef USE_MKL_VML
#include <mkl_vml_functions.h>
#endif

void FatVectorTransform::setBranchDependencies(const std::vector< std::vector<ForestNode*> >& aNodesByLevel)
{
	// Push only the branch id's (and compute num branches). The root is not pushed!
	mNumBranches = 0;
	mBranchByLevel.clear();
	std::vector< std::vector<ForestNode*> >::const_reverse_iterator inbl;
	for(inbl=aNodesByLevel.rbegin(); inbl != aNodesByLevel.rend(); ++inbl)
	{
		std::vector<unsigned int> v;
		std::vector<ForestNode*>::const_iterator ifn=inbl->begin();
		for(; ifn != inbl->end(); ++ifn)
		{
			v.push_back((*ifn)->mBranchId);
			++mNumBranches;
		}
		mBranchByLevel.push_back(v);
	}

	// Mark the first branch for nodes at the level below
	mFirstForLevel.assign(mNumBranches, false);
    ForestNode* curr_node = 0;
	for(inbl=aNodesByLevel.rbegin(); inbl != aNodesByLevel.rend(); ++inbl)
    {
		std::vector<ForestNode*>::const_iterator ifn=inbl->begin();
		for(; ifn != inbl->end(); ++ifn)
		{
			ForestNode* parent_node = (*ifn)->mParent;

			// If this is the first visit to the parent copy the result, otherwise do a element by element multiplication
			if(parent_node != curr_node)
			{
				curr_node = parent_node;
				mFirstForLevel[(*ifn)->mBranchId] = true;
			}
		}
	}

	// Discover the parent of each branch
	mParentNode.resize(mNumBranches);
	for(inbl=aNodesByLevel.rbegin(); inbl != aNodesByLevel.rend(); ++inbl)
	{
		std::vector<ForestNode*>::const_iterator ifn=inbl->begin();
		for(; ifn != inbl->end(); ++ifn)
		{
			mParentNode[(*ifn)->mBranchId] = (*ifn)->mParent->mBranchId+1; // The parent node is the node from which the branch originate
		}
	}
}


void FatVectorTransform::printCountGoodElements(void) const
{
	std::cout << std::endl;
	for(unsigned int b=0; b < mNumBranches; ++b)
	{
		size_t begin_idx = 0;
		for(; begin_idx < mNumSites; ++begin_idx)
		{
			int x = mNodeStatus[b*mNumSites+begin_idx];
			if(x == FatVectorTransform::SITE_EXISTS) break;
		}
		if(begin_idx == mNumSites)
		{
			std::ostringstream o;
			o << "No SITE_EXISTS in mNodePresent at branch: " << b << std::endl;
			throw FastCodeMLFatal(o);
		}

		size_t end_idx = mNumSites;
		for(; end_idx > begin_idx; --end_idx)
		{
			int x = mNodeStatus[b*mNumSites+end_idx-1];
			if(x == FatVectorTransform::SITE_EXISTS) break;
		}

		// Count the good elements
		unsigned int cnt = 0;
		for(size_t k=begin_idx; k < end_idx; ++k) if(mNodeStatus[b*mNumSites+k] == FatVectorTransform::SITE_EXISTS) ++cnt;

		std::cout << std::setw(2) << b << ": " << std::setw(4) << begin_idx << '-' << std::setw(4) << end_idx-1 << " (" << cnt << ")" << std::endl;
	}
}


void FatVectorTransform::printBranchVisitSequence(void) const
{
	std::cout << std::endl << "Branch at level" << std::endl;
	unsigned int level = 1;
	std::vector< std::vector<unsigned int> >::const_iterator inbl=mBranchByLevel.begin();
	const std::vector< std::vector<unsigned int> >::const_iterator end=mBranchByLevel.end();
	for(; inbl != end; ++inbl, ++level)
	{
		std::cout << level << ": ";

		std::vector<unsigned int>::const_iterator ifn=inbl->begin();
		for(; ifn != inbl->end(); ++ifn) std::cout << (*ifn) << ' ';

		std::cout << std::endl;
	}

	std::cout << std::endl << "Parent node for branch" << std::endl;
	for(unsigned int i=0; i < mNumBranches; ++i)
	{
		std::cout << std::setw(2) << i << " -> " << std::setw(2) << mParentNode[i] << std::endl;
	}
}



void FatVectorTransform::printNodeStatus(void) const
{
	std::cout << std::endl;
	for(unsigned int b=0; b < mNumBranches; ++b)
	{
		std::cout << "Branch " << b << std::endl;
		bool is_num = false;
		for(unsigned int k = 0; k < mNumSites; ++k)
		{
			int x = mNodeStatus[b*mNumSites+k];
			if(x == FatVectorTransform::SITE_NOT_EXISTS)  {std::cout << '-'; is_num = false;}
			else if(x == FatVectorTransform::SITE_EXISTS) {std::cout << 'x'; is_num = false;}
			else                                          {if(is_num) std::cout << ' '; std::cout << x; is_num = true;}
		}
		std::cout << std::endl << std::endl;
	}
}


void FatVectorTransform::compactMatrix(void)
{
	// For each branch
	unsigned int b;
	for(b=0; b < mNumBranches; ++b)
	{
		// Compute the index of the first valid site
		size_t begin_idx = 0;
		for(; begin_idx < mNumSites; ++begin_idx)
		{
			if(mNodeStatus[b*mNumSites+begin_idx] == FatVectorTransform::SITE_EXISTS) break;
		}
		if(begin_idx == mNumSites)
		{
			std::ostringstream o;
			o << "No SITE_EXISTS in mNodePresent at branch: " << b << std::endl;
			throw FastCodeMLFatal(o);
		}

		// Compute the last valid site (actually it points one after)
		size_t end_idx = mNumSites;
		for(; end_idx > begin_idx; --end_idx)
		{
			if(mNodeStatus[b*mNumSites+end_idx-1] == FatVectorTransform::SITE_EXISTS) break;
		}

		// Get the compaction moves
		VectorOfRanges cmds;
		for(int site_to=static_cast<int>(end_idx)-1; site_to >= static_cast<int>(begin_idx); --site_to)
		{
			// Select the first hole (from right)
			if(mNodeStatus[b*mNumSites+site_to] == FatVectorTransform::SITE_EXISTS) continue;

			// From left find the first valid entry
			unsigned int site_from = static_cast<unsigned int>(begin_idx);

			// Save the move command
			cmds.push_back(Range(site_from, site_to));

			// Update the left limit
			for(++begin_idx; begin_idx < (unsigned int)site_to; ++begin_idx)
			{
				// Select the first valid site (from left)
				if(mNodeStatus[b*mNumSites+begin_idx] == FatVectorTransform::SITE_EXISTS) break;
			}
		}
		mCopyCmds.push_back(cmds);

		// Save the new start index and count
		mLimits[b] = std::make_pair(begin_idx, end_idx-begin_idx);

		// Compute the reuse of another value moves
		VectorOfRangesNoCnt reuse;
		for(unsigned int k=0; k < mNumSites; ++k)
		{
			// Select a reuse pointer
			int x = mNodeStatus[b*mNumSites+k];
			if(x >= FatVectorTransform::SITE_FIRST_NUM)
			{
				reuse.push_back(RangeNoCnt(x, k));
			}
		}
		mReuseCmds.push_back(reuse);
	}

	// Remove the node status array no more needed
	//mNodeStatus.clear();
	std::vector<int> x;
	mNodeStatus.swap(x); // To really release memory

	// Try to combine contiguous ranges
	for(b=0; b < mNumBranches; ++b)
	{
		const size_t nc = mCopyCmds[b].size();
		if(nc < 2) continue;

		// Start with two valid
		for(size_t i=0; i < nc-1;)
		{
			if(mCopyCmds[b][i].from+1 == mCopyCmds[b][i+1].from && mCopyCmds[b][i].to == mCopyCmds[b][i+1].to+1)
			{
				// Try to extend the range to other with the same ordering
				size_t j=i+1;
				for(; j < nc-1; ++j)
				{
					if(mCopyCmds[b][j].from+1 != mCopyCmds[b][j+1].from || mCopyCmds[b][j].to != mCopyCmds[b][j+1].to+1) break;
				}

				// Update the command list
				// Example: (100, 10, 1) (101, 9, 1) --> (100, 9, 2) (101, 9, 0)
				mCopyCmds[b][i].cnt = static_cast<unsigned int>(j-i+1);
				mCopyCmds[b][i].to  = mCopyCmds[b][j].to;
				for(size_t k=i+1; k <= j; ++k) mCopyCmds[b][k].cnt = 0;

				i = j+1;
			}
			else
			{
				++i;
			}
		}
	}
}

void FatVectorTransform::printCommands(void) const
{
	for(unsigned int b=0; b < mNumBranches; ++b)
	{
		std::cout << std::endl << "*** Branch " << b << std::endl;

		VectorOfRanges::const_iterator icc=mCopyCmds[b].begin();
		const VectorOfRanges::const_iterator endc=mCopyCmds[b].end();
		for(; icc != endc; ++icc)
		{
			if(icc->cnt == 1)
				std::cout << "C " << std::setw(4) << icc->from << " - " << std::setw(4) << icc->to << std::endl;
			else if(icc->cnt > 1)
				std::cout << "C " << std::setw(4) << icc->from << " - " << std::setw(4) << icc->to << " (" << icc->cnt << ")" << std::endl;
		}

		VectorOfRangesNoCnt::const_iterator icr=mReuseCmds[b].begin();
		const VectorOfRangesNoCnt::const_iterator endr=mReuseCmds[b].end();
		for(; icr != endr; ++icr)
		{
			std::cout << "R " << std::setw(4) << icr->from << " - " << std::setw(4) << icr->to << std::endl;
		}

		std::cout << "L   from: " << mLimits[b].first << " cnt: " << mLimits[b].second << std::endl;
	}
}


void FatVectorTransform::preCompactLeaves(CacheAlignedDoubleVector& aProbs)
{
	// If the forest has not been reduced do nothing (this should not happens)
	if(mNoTransformations) return;

	// Find all the nodes that are parent of other nodes (ie. they are not leaves)
	std::set<unsigned int> non_leaves;
	non_leaves.insert(mParentNode.begin(), mParentNode.end());

	// Make a list of leaf nodes
	std::vector<unsigned int> leaves;
	for(unsigned int node=1; node <= mNumBranches; ++node)
	{
		// Check if the node is a leaf, if not skip it
		if(non_leaves.find(node) == non_leaves.end()) leaves.push_back(node);
	}

	// For all leaves and all sets
	const int len = static_cast<int>(leaves.size()*Nt);
#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(len, leaves, aProbs) schedule(guided)
#else
	#pragma omp parallel for default(shared) schedule(guided)
#endif
	for(int i=0; i < len; ++i)
	{
		const unsigned int node_idx = i / Nt;
		const unsigned int node     = leaves[node_idx];
		const unsigned int set_idx  = i-node_idx*Nt; // Was: i % Nt;
		const size_t       start    = VECTOR_SLOT*(mNumSites*Nt*node+set_idx*mNumSites);

		// Do all the copies as requested
		VectorOfRanges::const_iterator icc=mCopyCmds[node-1].begin();
		const VectorOfRanges::const_iterator end=mCopyCmds[node-1].end();
		for(; icc != end; ++icc)
		{
			if(icc->cnt == 1)
			{
				memcpy(&aProbs[start+VECTOR_SLOT*icc->to], &aProbs[start+VECTOR_SLOT*icc->from], N*sizeof(double));
			}
			else if(icc->cnt > 1)
			{
				memcpy(&aProbs[start+VECTOR_SLOT*icc->to], &aProbs[start+VECTOR_SLOT*icc->from], (VECTOR_SLOT*icc->cnt-(VECTOR_SLOT-N))*sizeof(double));
			}
		}
	}
}


void FatVectorTransform::postCompact(CacheAlignedDoubleVector& aStepResults, CacheAlignedDoubleVector& aProbs, unsigned int aLevel, unsigned int aNumSets)
{
	const int nsns = static_cast<int>(VECTOR_SLOT*mNumSites*aNumSets);
	if(mNoTransformations)
	{
		const size_t num_branch = mBranchByLevel[aLevel].size();
		for(size_t b=0; b < num_branch; ++b)
		{
			const unsigned int   my_branch = mBranchByLevel[aLevel][b];
			const unsigned int parent_node = mParentNode[my_branch];
			const unsigned int     my_node = my_branch + 1;

			if(mFirstForLevel[my_branch])
			{
				memcpy(&aProbs[VECTOR_SLOT*mNumSites*Nt*parent_node],
					   &aStepResults[VECTOR_SLOT*mNumSites*Nt*my_node],
					   (VECTOR_SLOT*mNumSites*aNumSets-(VECTOR_SLOT-N))*sizeof(double));
			}
			else
			{
#ifdef USE_MKL_VML
				const unsigned int start_parent = VECTOR_SLOT*mNumSites*Nt*parent_node;
				const unsigned int start_child  = VECTOR_SLOT*mNumSites*Nt*my_node;
				vdMul(nsns, &aProbs[start_parent], &aStepResults[start_child], &aProbs[start_parent]);
#else
#ifdef _MSC_VER
				#pragma omp parallel for default(none) shared(parent_node, my_node, aNumSets, aProbs, aStepResults, nsns) schedule(guided)
#else
				#pragma omp parallel for default(shared) schedule(guided)
#endif
                for(int i=0; i < nsns; ++i)
                {
                    aProbs[VECTOR_SLOT*mNumSites*Nt*parent_node+i] *= aStepResults[VECTOR_SLOT*mNumSites*Nt*my_node+i];
                }
#endif
			}
		}
	}
	else
	{
		// For all the branches just processed
		const size_t num_branch = mBranchByLevel[aLevel].size();
		for(size_t b=0; b < num_branch; ++b)
		{
			const unsigned int branch      = mBranchByLevel[aLevel][b];
			const unsigned int node        = branch + 1;
			const unsigned int parent_node = mParentNode[branch];

			// Reverse all copies (copy back the values copied in the previous step to fill holes)
			VectorOfRanges::const_iterator icc=mCopyCmds[branch].begin();
			const VectorOfRanges::const_iterator end=mCopyCmds[branch].end();
			for(; icc != end; ++icc)
			{
				// Make local copies to increase locality
				unsigned int cnt = icc->cnt;

				if(cnt == 1)
				{
					const size_t from_idx = VECTOR_SLOT*(mNumSites*Nt*node+icc->from);
					const size_t to_idx   = VECTOR_SLOT*(mNumSites*Nt*node+icc->to);

					for(unsigned int set_idx=0; set_idx < aNumSets; ++set_idx)
					{
						memcpy(&aStepResults[from_idx+set_idx*mNumSites*VECTOR_SLOT],
							   &aStepResults[to_idx+set_idx*mNumSites*VECTOR_SLOT],
							   N*sizeof(double));
					}
				}
				else if(cnt > 1)
				{
					const size_t from_idx = VECTOR_SLOT*(mNumSites*Nt*node+icc->from);
					const size_t to_idx   = VECTOR_SLOT*(mNumSites*Nt*node+icc->to);

					for(unsigned int set_idx=0; set_idx < aNumSets; ++set_idx)
					{
						memcpy(&aStepResults[from_idx+set_idx*mNumSites*VECTOR_SLOT],
							   &aStepResults[to_idx+set_idx*mNumSites*VECTOR_SLOT],
							   (VECTOR_SLOT*cnt-(VECTOR_SLOT-N))*sizeof(double));
					}
				}
			}

			// Reuse values 
			VectorOfRangesNoCnt::const_iterator icr=mReuseCmds[branch].begin();
			const VectorOfRangesNoCnt::const_iterator endr=mReuseCmds[branch].end();
			for(; icr != endr; ++icr)
			{
				// Make local copies to increase locality
				const size_t from_idx = VECTOR_SLOT*(mNumSites*Nt*node+icr->from);
				const size_t to_idx   = VECTOR_SLOT*(mNumSites*Nt*node+icr->to);

				for(unsigned int set_idx=0; set_idx < aNumSets; ++set_idx)
				{
					memcpy(&aStepResults[to_idx+set_idx*mNumSites*VECTOR_SLOT],
						   &aStepResults[from_idx+set_idx*mNumSites*VECTOR_SLOT],
						   N*sizeof(double));
				}
			}

			if(mFirstForLevel[branch])
			{
				memcpy(&aProbs[VECTOR_SLOT*mNumSites*Nt*parent_node], &aStepResults[VECTOR_SLOT*mNumSites*Nt*node], VECTOR_SLOT*mNumSites*aNumSets*sizeof(double));
			}
			else
			{
#ifdef USE_MKL_VML
				const unsigned int start_parent = VECTOR_SLOT*mNumSites*Nt*parent_node;
				const unsigned int start_child  = VECTOR_SLOT*mNumSites*Nt*node;
				vdMul(nsns, &aProbs[start_parent], &aStepResults[start_child], &aProbs[start_parent]);
#else
#ifdef _MSC_VER
				#pragma omp parallel for default(none) shared(parent_node, node, aNumSets, aProbs, aStepResults, nsns) schedule(guided)
#else
				#pragma omp parallel for default(shared) schedule(guided)
#endif
                for(int i=0; i < nsns; ++i)
                {
                    aProbs[VECTOR_SLOT*mNumSites*Nt*parent_node+i] *= aStepResults[VECTOR_SLOT*mNumSites*Nt*node+i];
                }
#endif
			}

			// Copy for the next branch (if this branch does not lead to the root)
			if(parent_node)
			{
				// Do all the copies as requested
				for(icc=mCopyCmds[parent_node-1].begin(); icc != mCopyCmds[parent_node-1].end(); ++icc)
				{
					// Make local copies to increase locality
					unsigned int cnt = icc->cnt;

					if(cnt == 1)
					{
						const size_t from_idx = VECTOR_SLOT*(mNumSites*Nt*parent_node+icc->from);
						const size_t to_idx   = VECTOR_SLOT*(mNumSites*Nt*parent_node+icc->to);

						for(unsigned int set_idx=0; set_idx < aNumSets; ++set_idx)
						{
							memcpy(&aProbs[to_idx+set_idx*mNumSites*VECTOR_SLOT],
								   &aProbs[from_idx+set_idx*mNumSites*VECTOR_SLOT],
								   N*sizeof(double));
						}
					}
					if(cnt > 1)
					{
						const size_t from_idx = VECTOR_SLOT*(mNumSites*Nt*parent_node+icc->from);
						const size_t to_idx   = VECTOR_SLOT*(mNumSites*Nt*parent_node+icc->to);

						for(unsigned int set_idx=0; set_idx < aNumSets; ++set_idx)
						{
							memcpy(&aProbs[to_idx+set_idx*mNumSites*VECTOR_SLOT],
								   &aProbs[from_idx+set_idx*mNumSites*VECTOR_SLOT],
								   (VECTOR_SLOT*cnt-(VECTOR_SLOT-N))*sizeof(double));
						}
					}
				}
			}
		}
	}
}

#endif
