/// @mainpage FastCodeML
///
/// @section intro_sect Introduction
/// 
/// FastCodeML is a rewrite of CodeML based directly on the pseudocode document.
/// It incorporates various parallelization strategies to be able to exploit modern HPC machines architecture.
/// For this reason there are various parts of the code that can be selected at compile time or run time to experiment with various, possible solutions.
///
/// @section contacts_sect Contacts
///
/// Contact us if you want more information on the project, want to collaborate or suggest new ideas.
/// 
///- Ing. <a href="mailto:mvalle@cscs.ch">Mario Valle</a> - Swiss National Supercomputing Centre (CSCS) - Switzerland
///- The HP2C <a href="mailto:selectome@hp2c.ch">Selectome</a> Project Group - Mainly based in University of Lausanne - Switzerland
///

#include <iostream>
#include <iomanip>
#include <limits>
#include "CmdLine.h"
#include "Newick.h"
#include "Phylip.h"
#include "BayesTest.h"
#include "Forest.h"
#include "Exceptions.h"
#include "BranchSiteModel.h"
#include "ParseParameters.h"
#include "VerbosityLevels.h"
#include "WriteResults.h"

#ifndef VTRACE
#ifdef _OPENMP
#include <omp.h>
#endif
#endif
#ifdef USE_MKL_VML
#include <mkl_vml.h>
#endif
#include "Timer.h"
#ifdef USE_MPI
#include "HighLevelCoordinator.h"
#endif

/// Main program for FastCodeML.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-12-22 (initial version)
///     @version 1.0
///
///	@param[in] aRgc Number of command line parameters
/// @param[in] aRgv Command line parameters
///
int main(int aRgc, char **aRgv)
{
	try
	{
#ifdef USE_MKL_VML
	// If used, intitialize the MKL VML library
	vmlSetMode(VML_HA|VML_DOUBLE_CONSISTENT);
#endif

#ifdef USE_MPI
	// Start the high level parallel executor (based on MPI)
	HighLevelCoordinator hlc(&aRgc, &aRgv);
#endif

	// Parse the command line
	CmdLine cmd;
	cmd.parseCmdLine(aRgc, aRgv);

	// Adjust and report the number of threads that will be used
#ifdef _OPENMP
	int num_threads = omp_get_max_threads();
	if(num_threads < 2 || cmd.mForceSerial)
	{
		cmd.mForceSerial = true;
		num_threads = 1;
		omp_set_num_threads(1);
	}
#else
	cmd.mForceSerial = true;
	int num_threads = 1;
#endif

#ifdef USE_MPI
	// Shutdown messages from all MPI processes except the master
	if(!hlc.isMaster()) cmd.mVerboseLevel = VERBOSE_NONE;
#endif

	// Write out command line parameters (if not quiet i.e. if verbose level > 0)
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT)
	{
													std::cout << std::endl;
													std::cout << "Tree file:      " << cmd.mTreeFile << std::endl;
													std::cout << "Gene file:      " << cmd.mGeneFile << std::endl;
													std::cout << "Verbose level:  " << cmd.mVerboseLevel << " (" << decodeVerboseLevel(cmd.mVerboseLevel) << ')' << std::endl;
		if(cmd.mSeed)								std::cout << "Seed:           " << cmd.mSeed << std::endl;
		if(cmd.mBranchFromFile)						std::cout << "Branch:         From tree file" << std::endl;
		else if(cmd.mBranchStart != UINT_MAX && cmd.mBranchStart == cmd.mBranchEnd)
			                                        std::cout << "Branch:         " << cmd.mBranchStart << std::endl;
		else if(cmd.mBranchStart != UINT_MAX && cmd.mBranchEnd == UINT_MAX)
			                                        std::cout << "Branches:       " << cmd.mBranchStart << "-end" << std::endl;
		else if(cmd.mBranchStart != UINT_MAX && cmd.mBranchEnd != UINT_MAX)
													std::cout << "Branches:       " << cmd.mBranchStart << '-' << cmd.mBranchEnd << std::endl;
		if(!cmd.mStopIfNotLRT)						std::cout << "H0 pre stop:    No" << std::endl;
		if(cmd.mIgnoreFreq)							std::cout << "Codon freq.:    Ignore" << std::endl;
		if(cmd.mDoNotReduceForest)					std::cout << "Reduce forest:  Do not reduce" << std::endl;
		else										std::cout << "Reduce forest:  Aggressive" << std::endl;
		if(cmd.mInitH0fromH1)						std::cout << "Starting val.:  From H1" << std::endl;
		else if(cmd.mInitFromParams && cmd.mBranchLengthsFromFile)
													std::cout << "Starting val.:  Times from tree file and params from const (see below)" << std::endl;
		else if(cmd.mInitFromParams)				std::cout << "Starting val.:  Params from const (see below)" << std::endl;
		else if(cmd.mBranchLengthsFromFile)			std::cout << "Starting val.:  Times from tree file" << std::endl;
		if(cmd.mNoMaximization)						std::cout << "Maximization:   No" << std::endl;
		if(cmd.mTrace)								std::cout << "Trace:          On" << std::endl;
		if(cmd.mCleanData)							std::cout << "Clean data:     On" << std::endl;
		else										std::cout << "Clean data:     Off" << std::endl;
		if(cmd.mGraphFile)							std::cout << "Graph file:     " << cmd.mGraphFile << std::endl;
		if(cmd.mGraphFile && cmd.mExportComputedTimes != UINT_MAX)
													std::cout << "Graph times:    From H" << cmd.mExportComputedTimes << std::endl;
		if(!cmd.mNoMaximization)					std::cout << "Optimizer:      " << cmd.mOptimizationAlgo << std::endl;
		if(cmd.mMaxIterations != MAX_ITERATIONS)	std::cout << "Max iterations: " << cmd.mMaxIterations << std::endl;
		if(cmd.mDeltaValueForGradient > 0.0)		std::cout << "Delta value:    " << cmd.mDeltaValueForGradient << std::endl;
													std::cout << "Relative error: " << cmd.mRelativeError << std::endl;
		if(cmd.mResultsFile)						std::cout << "Results file:   " << cmd.mResultsFile << std::endl;

#ifdef _OPENMP
		if(num_threads > 1)
		{
													std::cout << "Num. threads:   " << num_threads << std::endl
		                                                      << "Num. cores:     " << omp_get_num_procs() << std::endl;
		}
		else
#endif
		{
													std::cout << "Num. threads:   1 serial" << std::endl
		                                                      << "Num. cores:     1"  << std::endl;
		}
#ifdef USE_MPI
		if(hlc.numJobs() > 2)						std::cout << "Num. MPI proc:  1 (master) + " << hlc.numJobs()-1 << " (workers)" << std::endl;
		else										std::cout << "Num. MPI proc:  Insufficient, single task execution" << std::endl;
#endif
													std::cout << "Compiled with:  ";
#ifdef _OPENMP
													std::cout << "USE_OPENMP ";
#endif
#ifdef USE_MPI
													std::cout << "USE_MPI ";
#endif
#ifdef USE_CPV_SCALING
													std::cout << "USE_CPV_SCALING ";
#endif
#ifdef NEW_LIKELIHOOD
													std::cout << "NEW_LIKELIHOOD ";
#endif
#ifdef NON_RECURSIVE_VISIT
													std::cout << "NON_RECURSIVE_VISIT ";
#endif
#ifdef USE_DAG
													std::cout << "USE_DAG ";
#endif
#ifdef USE_ORIGINAL_PROPORTIONS
													std::cout << "USE_ORIGINAL_PROPORTIONS ";
#endif
#ifdef USE_LAPACK
													std::cout << "USE_LAPACK ";
#endif
#ifdef USE_MKL_VML
													std::cout << "USE_MKL_VML";
#endif
													std::cout << std::endl << std::endl;
													if(cmd.mInitFromParams)
													{
														std::cout << "Param initial values:" << std::endl << std::endl
																  << ParseParameters::getInstance();
													}
	}

	// Initialize the random number generator (0 means it is not set on the command line)
#ifdef USE_MPI
	// Insure that each MPI process starts with a different seed
	if(cmd.mSeed == 0) cmd.mSeed = static_cast<unsigned int>(time(NULL)) + static_cast<unsigned int>(hlc.getRank()) * 1000;
#else
	if(cmd.mSeed == 0) cmd.mSeed = static_cast<unsigned int>(time(NULL));
#endif
	srand(cmd.mSeed);

	// Verify the optimizer algorithm selected on the command line
	if(!cmd.mNoMaximization) BranchSiteModel::verifyOptimizerAlgo(cmd.mOptimizationAlgo);

	// Start a timer (to measure serial part over parallel one)
	Timer timer;
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) timer.start();

	// Create the forest
	Forest forest(cmd.mVerboseLevel);

	// Enclose file loading into a block so temporary structures could be deleted when no more needed
	{
	// Load the multiple sequence alignment (MSA)
	Phylip msa(cmd.mVerboseLevel);
	msa.readFile(cmd.mGeneFile, cmd.mCleanData);

	// Load the phylogenetic tree
	Newick tree(cmd.mVerboseLevel);
	tree.readFile(cmd.mTreeFile);

	// Check coherence between the two files
	msa.checkNameCoherence(tree.getSpecies());

	// Check root
	tree.checkRootBranches();

	// If times from file then check for null branch lengths for any leaf
	if(cmd.mBranchLengthsFromFile)
	{
		int zero_on_leaf_cnt = 0;
		int zero_on_int_cnt  = 0;
		tree.countNullBranchLengths(zero_on_leaf_cnt, zero_on_int_cnt);

		if(zero_on_leaf_cnt > 0 || zero_on_int_cnt > 0)
		{
			std::cout << "Found null or missing branch length in tree file. On leaves: " << zero_on_leaf_cnt << "  on internal branches: " << zero_on_int_cnt << std::endl;
		}

		if(zero_on_leaf_cnt > 0)
		{
			throw FastCodeMLFatal("Null or missing branch length in tree file");
		}
	}

	// Print the tree with the numbering of internal branches
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) tree.printTreeAnnotated(std::cout);

	// Load the forest
	forest.loadTreeAndGenes(tree, msa, cmd.mIgnoreFreq ? CodonFrequencies::CODON_FREQ_MODEL_UNIF : CodonFrequencies::CODON_FREQ_MODEL_F3X4);
	}

	// Reduce the forest merging common subtrees. Add also more reduction, then clean the no more useful data.
	if(!cmd.mDoNotReduceForest)
	{
		//bool sts = forest.reduceSubtrees(cmd.mNumReductionBlocks);
		forest.reduceSubtrees();

#ifndef NEW_LIKELIHOOD
		forest.addAggressiveReduction();
#endif
		forest.cleanReductionWorkingData();		
#ifdef NEW_LIKELIHOOD
		forest.prepareNewReduction();
#endif
	}
#ifdef NEW_LIKELIHOOD
	else
	{
		forest.prepareNewReductionNoReuse();
	}
#endif

#ifdef NON_RECURSIVE_VISIT
	// Prepare the pointers to visit the trees without recursion
	forest.prepareNonRecursiveVisit();
#endif

	// Subdivide the trees in groups based on dependencies
	//forest.prepareDependencies(cmd.mForceSerial || cmd.mDoNotReduceForest);

#ifdef USE_DAG
	// Load the forest into a DAG
	forest.loadForestIntoDAG(Nt);
#endif

	// Get the time needed by data preprocessing
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) {timer.stop(); std::cout << std::endl << "TIMER (preprocessing) ncores: " << std::setw(2) << num_threads << " time: " << timer.get() << std::endl;}

	// Print few statistics
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) std::cout << forest;

#ifdef USE_MPI
	// Distribute the work. If run under MPI then finish, else return to the standard execution flow
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) timer.start();
	bool has_run_under_MPI = hlc.startWork(forest, cmd);

	// If executed under MPI report the time spent, otherwise stop the timer so it can be restarted around the serial execution
	if(has_run_under_MPI)
	{
		if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) {timer.stop(); std::cout << std::endl << "TIMER (processing) ncores: " << std::setw(2) << num_threads*(hlc.numJobs()-1)+1 << " time: " << timer.get() << std::endl;}
		return 0;
	}
	else
	{
		timer.stop();
	}
#endif

	// Initialize the output results file (if the argument is null, no file is created)
	WriteResults output_results(cmd.mResultsFile);

	// Compute the range of branches to mark as foreground
	size_t branch_start, branch_end;
	forest.getBranchRange(cmd, branch_start, branch_end);

	// Start timing parallel part
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) timer.start();

	// Initialize the models
	BranchSiteModelNullHyp h0(forest, cmd);
	BranchSiteModelAltHyp  h1(forest, cmd);

	// Initialize the test
	BayesTest beb(forest, cmd.mVerboseLevel, cmd.mDoNotReduceForest);

	// For all requested internal branches
	for(size_t fg_branch=branch_start; fg_branch <= branch_end; ++fg_branch)
	{
		if(cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) std::cout << std::endl << "Doing branch " << fg_branch << std::endl;

		// Compute the alternate model maximum loglikelihood
		double lnl1 = 0.;
		if(cmd.mComputeHypothesis != 0)
		{
			if(cmd.mInitFromParams)			h1.initFromParams();
			if(cmd.mBranchLengthsFromFile)	h1.initFromTree();

			lnl1 = h1(fg_branch);

			// Save the value for formatted output
			output_results.saveLnL(fg_branch, lnl1, 1);
		}

		// Compute the null model maximum loglikelihood
		double lnl0 = 0.;
		if(cmd.mComputeHypothesis != 1)
		{
			if(cmd.mInitH0fromH1)				h0.initFromResult(h1.getVariables());
			else
			{
				if(cmd.mInitFromParams)			h0.initFromParams();
				if(cmd.mBranchLengthsFromFile)	h0.initFromTree();
			}

			lnl0 = h0(fg_branch, cmd.mStopIfNotLRT && cmd.mComputeHypothesis != 0, lnl1-THRESHOLD_FOR_LRT);

			// Save the value for formatted output (only if has not be forced to stop)
			if(lnl0 < DBL_MAX) output_results.saveLnL(fg_branch, lnl0, 0);
		}

		if(cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS)
		{
			std::cout << std::endl;
			if(cmd.mComputeHypothesis != 1)
			{
				std::cout << "LnL0: ";
				if(lnl0 == std::numeric_limits<double>::infinity())
					std::cout << "**Invalid result**";
				else if(lnl0 < DBL_MAX)
					std::cout << std::setprecision(15) << std::fixed << lnl0;
				else
					std::cout << "(Doesn't pass LRT, skipping)";
				std::cout << " Function calls: " << h0.getNumEvaluations() << "   ";
				std::cout << std::endl << std::endl;
				if(lnl0 != std::numeric_limits<double>::infinity()) h0.printFinalVars(std::cout);
				std::cout << std::endl;
			}
			if(cmd.mComputeHypothesis != 0)
			{
				std::cout << "LnL1: ";
				if(lnl1 == std::numeric_limits<double>::infinity())
					std::cout << "**Invalid result**";
				else
					std::cout << std::setprecision(15) << std::fixed << lnl1;
				std::cout << " Function calls: " << h1.getNumEvaluations();
				std::cout << std::endl << std::endl;
				if(lnl1 != std::numeric_limits<double>::infinity()) h1.printFinalVars(std::cout);
				std::cout << std::endl;
			}
			if(cmd.mComputeHypothesis > 1)
			{
				if(lnl0 == std::numeric_limits<double>::infinity() || lnl1 == std::numeric_limits<double>::infinity())
					std::cout << "LRT: **Invalid result**";
				else if(lnl0 < DBL_MAX)
					std::cout << "LRT: " << std::setprecision(15) << std::fixed << lnl1 - lnl0 << "  (threshold: " << std::setprecision(15) << std::fixed << THRESHOLD_FOR_LRT << ')';
				else
					std::cout << "LRT: < " << std::setprecision(15) << std::fixed << THRESHOLD_FOR_LRT;
				std::cout << std::endl;
			}
		}

		// If requested set the time in the forest and export to a graph visualization tool
		if(cmd.mGraphFile)
		{
			switch(cmd.mExportComputedTimes)
			{
			case 0:
				h0.saveComputedTimes();
				break;

			case 1:
				h1.saveComputedTimes();
				break;

			default:
				break;
			}

			// Use the forest export class
			ForestExport fe(forest);
			fe.exportForest(cmd.mGraphFile, fg_branch);
		}

		// If the two hypothesis are computed, H0 has not been stopped and the run passes the LRT, then compute the BEB
		if(cmd.mComputeHypothesis > 1 && lnl0 < DBL_MAX && BranchSiteModel::performLRT(lnl0, lnl1))
		{
			// Get the scale values from the latest optimized h1.
			std::vector<double> scales(2);
			h1.getScales(scales);

			// Run the BEB test
			beb.computeBEB(h1.getVariables(), fg_branch, scales);

			// Output the sites under positive selection (if any)
			if(cmd.mVerboseLevel >= VERBOSE_ONLY_RESULTS) beb.printPositiveSelSites(fg_branch);

			// Get the sites under positive selection for printing in the results file (if defined)
			if(output_results.isWriteResultsEnabled())
			{
				std::vector<unsigned int> positive_sel_sites;
				std::vector<double> positive_sel_sites_prob;
				beb.extractPositiveSelSites(positive_sel_sites, positive_sel_sites_prob);
				output_results.savePositiveSelSites(fg_branch, positive_sel_sites, positive_sel_sites_prob);
			}
		}
	}

	// Get the time needed by the parallel part
	if(cmd.mVerboseLevel >= VERBOSE_INFO_OUTPUT) {timer.stop(); std::cout << std::endl << "TIMER (processing) ncores: " << std::setw(2) << num_threads << " time: " << timer.get() << std::endl;}

	// Output the results
	output_results.outputResults();

	////////////////////////////////////////////////////////////////////
	// Catch all exceptions
	}
	catch(const FastCodeMLSuccess&)
	{
		return 0;
	}
	catch(const FastCodeMLFatal& e)
	{
		// If a message associated (i.e. no empty string), display it
		if(e.what()[0]) std::cout << std::endl << e.what() << std::endl << std::endl;
		return 1;
	}
	catch(const FastCodeMLMemoryError& e)
	{
		std::cout << std::endl << e.what() << std::endl << std::endl;
		return 1;
	}
	catch(const std::bad_alloc& e)
	{
		std::cout << std::endl << e.what() << std::endl << std::endl;
		return 1;
	}
	catch(...)
	{
		std::cout << std::endl << "Default exception caught." << std::endl << std::endl;
		return 1;
	}

	return 0;
}

/// @page cppstd_page C++ Coding Standard
/// Here are collected few rules for coding this project.
///
/// @section cnames_sect Class names
/// Class names are CamelCase with first letter uppercase.
///
/// Ex: %PhyloTree
///
/// @section cmeth_sect Class methods
/// Class methods names are CamelCase with the first letter lowercase.
/// Only very common and specific names should be all lowercase, like read, clean, size.
///
/// Ex: testFillQ
///
/// @section cdatamemb_sect Class data members
/// Class member variables names start with 'm' followed by CamelCase name.
///
/// Ex: mFgBranch
///
/// @section carg_sect Function arguments
/// Function arguments names start with 'a' followed by CamelCase name.
///
/// Ex: aFgBranch
///
/// @section const_sect Constants and enumeration
/// Constants and enumerations are all uppercase with words separated by '_'.
/// The first letters specify the kind of constant (like: STS_ for status, OPT_ for option value).
///
/// Ex: STS_CANT_OPEN
///
/// @section stack_sect Temporary variables
/// All the other variables are all lower case with parts separated by '_'.
///
/// Ex: branch_list
///
/// @section misc_sect Miscellaneous rules
/// In case of error main should return 1.
///
/// Array sizes and corresponding indexes should be size_t. The remaining counters should be unsigned int. 
///
/// The null pointer should be written as NULL, not 0 to make clear its purpose.
///

/**
@page cmd_page Command Line Switches
Here is a quick list of the valid command line switches for FastCodeML.

The input `tree_file` is in %Newick format with the file containing only one tree. The `alignment_file` instead is in %Phylip format.

@verbatim

Usage:
    FastCodeML [options] tree_file alignment_file

-d  --debug  -v  --verbose (required argument)
        Verbosity level (0: none; 1: results only; 2: normal info; 3: MPI trace; 4: more debug) (default: 1)

-q  --quiet (no argument)
        No messages except results

-?  -h  --help (no argument)
        This help

-s  --seed (required argument)
        Random number generator seed (0 < seed < 1000000000)

-b  --branch (required argument)
        Do only this branch as foreground branch (count from 0)

-bs  --branch-start (required argument)
        Start computing from this branch as foreground one (count from 0) (default: first one)

-be  --branch-end (required argument)
        End computing at this branch as foreground one (count from 0) (default: last one)

-i  --ignore-freq (no argument)
        Ignore computed codon frequency and set all of them to 1/61

-e  --export (required argument)
        Export forest in GML format (if %03d or @03d is present, one is created for each fg branch)

-nr  --no-reduce (no argument)
        Do not reduce forest by merging common subtrees

-l  --lengths-from-file  --times-from-file (no argument)
        Initial branch lengths from tree file

-o  --initial-step (no argument)
        Only the initial step is performed (no maximization)

-c  --export-comp-times (required argument)
        Export the computed times from H0 if 0, H1 if 1, otherwise the one read in the phylo tree

-r  --trace (no argument)
        Trace the maximization run

-np  --no-parallel (no argument)
        Don't use parallel execution

-bf  --branch-from-file (no argument)
        Do only the branch marked in the file as foreground branch

-hy  --only-hyp (required argument)
        Compute only H0 if 0, H1 if 1

-i1  --init-from-h1 (no argument)
        Start H0 optimization from H1 results

-m  --maximizer (required argument)
        Optimizer algorithm (0:LBFGS, 1:VAR1, 2:VAR2, 3:SLSQP, 11:BOBYQA, 22:FromCodeML, 99:MLSL_LDS) (default: 0)

-sd  --small-diff (required argument)
        Delta used in gradient computation (default: 1.49e-8)

-p  --init-param (required argument)
        Pass initialization parameter in the form: P=value (P: w0, k, p0, p1, w2)

-ic  --init-default (no argument)
        Start from default parameter values and times from tree file

-x  --extra-debug (required argument)
        Extra debug parameter (zero disables it)

-re  --relative-error (required argument)
        Relative error where to stop maximization (default: 1e-3)

-ou  --output (required argument)
        Write results formatted to this file

-cl  --clean-data (no argument)
        Remove ambiguous or missing sites from the MSA (default: no)

-ps  --no-pre-stop (no argument)
        Don't stop H0 maximization even if it cannot satisfy LRT (default: stop)

-mi  --max-iterations (required argument)
        Maximum number of iterations for the maximizer (default: 10000)

@endverbatim
*/

/// @page vampir_page Using Vampir for profiling
/// On Linux we use VampirTrace to collect profile data and Vampir to display the results (http://www.vampir.eu/).
///
/// Before invoking CMAKE define CXX=vtCC
///
/// Define CMAKE_BUILD_TYPE as: RelWithDebInfo
///
/// Run CMAKE and configure.
///
/// Then define CMAKE_CXX_FLAGS as: -vt:mt -vt:noopari
/// If you want to trace also the OpenMP calls then change it to: -vt:mt -vt:preprocess -DVTRACE
///
/// Then proceed as usual to build the executable.
///
/// Before running the executable, define the following environment variables:
///
///     export VT_BUFFER_SIZE=512M
///     export VT_MAX_FLUSHES=0
///     export VT_SYNCH_FLUSH=yes
///     export VT_GPUTRACE=no
///     export VT_UNIFY=no
///
/// Due to a VampirTrace bug, at the end of the execution, run the vtunify executable by itself.
///
/// Now you can analyze the results by running vampir on the *.otf file generated.
///
