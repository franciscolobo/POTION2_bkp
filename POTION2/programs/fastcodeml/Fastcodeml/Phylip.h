
#ifndef PHYLIP_H
#define PHYLIP_H

#include <vector>
#include <string>
#include "Genes.h"

/// The genes of the set of species under analysis.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-02-15 (initial version)
///     @version 1.0
///
///
class Phylip : public Genes
{
public:
	/// Constructor
	///
	/// @param[in] aVerboseLevel The verbosity level
	///
	explicit Phylip(unsigned int aVerboseLevel=0) : Genes(aVerboseLevel) {}

private:
	/// Load the gene file in Phylip format.
	///
	/// @param[in] aFilename The filename containing the genes under analysis
	/// @param[out] aSpecies The list of species read
	/// @param[out] aSequences Array of genes, one per specie
	///
	/// @exception FastCodeMLFatal On various error conditions
	///
	virtual void loadData(const char* aFilename, std::vector<std::string>& aSpecies, std::vector<std::string>& aSequences);
};

#endif

