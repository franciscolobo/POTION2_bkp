
#ifndef GENES_H
#define GENES_H

#include <string>
#include <vector>
#include <map>

/// The genes of the set of species under analysis.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-12-22 (initial version)
///     @version 1.0
///
class Genes
{
protected:
	/// Constructor
	///
	/// @param[in] aVerboseLevel The verbosity level
	///
	explicit Genes(unsigned int aVerboseLevel=0) : mVerboseLevel(aVerboseLevel), mOriginalNumSites(0)
	{
		initFullCodonMap();
	}

	/// Destructor
	///
	virtual ~Genes();

public:
	/// Load the gene file.
	///
	/// @param[in] aFilename The filename containing the genes under analysis
	/// @param[in] aCleanData If true remove missing or ambiguous sites
	///
	/// @exception FastCodeMLFatal If cannot open file and other problems
	///
	void readFile(const char* aFilename, bool aCleanData);

	/// Return the number of valid sites loaded.
	///
	/// @return The number of sites
	///
	size_t getNumSites(void) const {return mSiteMultiplicity.size();}

	/// Return the site multiplicity.
	///
	/// @return Reference to the site multiplicity array
	///
	const std::vector<unsigned int>& getSiteMultiplicity(void) const {return mSiteMultiplicity;}

	/// Return the codon index for the one identified by the specie label and the site.
	/// It record internally the decoded codon to be used by setLeaveProb() and updateCodonCount()
	///
	/// @param[in] aSpecie The specie label.
	/// @param[in] aSite The site index.
	///
	/// @return The codon index or -1 in case of error or invalid codon at the given position.
	///
	long long getCodonIdx(std::string aSpecie, size_t aSite) const;

	/// Set the correct positions in the leave probability vector to 1/num_positions.
	///
	/// @param[out] aLeaveProbVect The leave probability vector to be set.
	///
	/// @exception FastCodeMLFatal If saved codon is invalid.
	///
	void setLeaveProb(double* aLeaveProbVect) const;

	/// Save codons in a given array for later count.
	/// For each codon the aCodon array contains one more vector that starts with the aSiteMultiplicity value followed by the codon positions
	/// that is, one position for non-ambiguous codons, 2 to 61 for ambiguous codons.
	///
	/// @param[in,out] aCodons Codon positions and multiplicity.
	/// @param[in] aSiteMultiplicity The multiplicity of the given site.
	///
	/// @exception FastCodeMLFatal If saved codon is invalid.
	///
	void saveCodonsForCount(std::vector<std::vector<unsigned int> >& aCodons, unsigned int aSiteMultiplicity) const;

	/// Check coherence between tree and genes.
	///
	/// @param[in] aNames The phylogenetic tree species names
	///
	/// @exception FastCodeMLFatal Throw exception if the species do not match
	///
	void checkNameCoherence(const std::vector<std::string>& aNames) const;

	/// Access the map that convert from the site number to the original MSA site number.
	///
	/// @return The map from reduced site number to list of corresponding original sites.
	///
	const std::multimap<size_t, size_t>& getSitesMappingToOriginal(void) const {return mSitesMappingToOriginal;}

	/// Get the number of sites in the loaded MSA.
	///
	/// @return The number of sites in the loaded gene file.
	///
	size_t getOriginalNumSites(void) const {return mOriginalNumSites;}

	/// Get the number of species in the loaded MSA.
	///
	/// @return The number of species in the loaded gene file.
	///
	size_t getNumSpecies(void) const {return mDnaSpecies.size();}


private:
	/// Convert given codon in the corresponding codons positions after expanding ambiguous characters.
	/// The valid letter are TGKCYSBAWRDMHVN and U that is mapped to T and - that is mapped to N.
	/// The full list can be found on http://en.wikipedia.org/wiki/Nucleic_acid_notation
	///
	/// @param[in] aCodon The three letters for the codon (no need to be zero terminated)
	///
	/// @return The list of expanded positions or an empty list on error.
	///
	const std::vector<int>& getPositions(const char* aCodon) const;

	/// Test if the three letters of the argument represent a valid codon.
	///
	/// @param[in] aCodon String of three letters representing the codon (not null terminated)
	/// @param[in] aRemoveAmbiguous If true only TCAG codons are valid, else also ambiguous one are valid
	///
	/// @return True if codon is a valid codon
	///
	bool validCodon(const char* aCodon, bool aRemoveAmbiguous) const;

	/// Load the gene file.
	/// This routine is redefined in every derived class to load a specific format.
	///
	/// @param[in] aFilename The filename containing the genes under analysis
	/// @param[out] aSpecies The list of species read
	/// @param[out] aSequences Array of genes, one per specie
	///
	/// @exception FastCodeMLFatalNoMsg On various error conditions
	///
	virtual void loadData(const char* aFilename, std::vector<std::string>& aSpecies, std::vector<std::string>& aSequences) =0;

	/// Initialize the valid codon map.
	/// This routine fills mMapCodonToPosition.
	///
	void initFullCodonMap(void);

	/// Compare two codons.
	///
	/// @param[in] aCodon1 First codon to be compared (three characters, no need for zero termination)
	/// @param[in] aCodon2 Second codon to be compared (three characters, no need for zero termination)
	///
	/// @return True if the codons are valid, equal or if ambiguous expand to the same set of positions.
	///
	bool compareCodons(const char* aCodon1, const char* aCodon2) const;

protected:
	unsigned int								mVerboseLevel;				///< The verbosity level as set in the constructor

private:
	std::vector<std::string>					mDnaSpecies;				///< The list of species labels
	std::vector<std::string>					mDnaGene;					///< The gene DNA basis strings
	std::vector<unsigned int>					mSiteMultiplicity;			///< Site multiplicity (sites with multiplicity of zero has been removed from the site list)
	std::vector<unsigned int>					mMapSiteToDnaGene;			///< Map the site number to the position in mDnaGene
	std::map<std::string, unsigned int>			mMapSpecieToDnaGene;		///< Map specie name to position in the gene list mDnaGene
	std::multimap<size_t, size_t>				mSitesMappingToOriginal;	///< Map reduced site num. to list of corresponding original sites
	size_t										mOriginalNumSites;			///< Original number of sites (before cleaning)

	std::map<std::string, std::vector<int> >	mMapCodonToPosition;		///< Map codons (including ambiguous ones) to positions on the CPV
	std::vector<int>							mEmptyVector;				///< Empty vector to be returned if no position available
	mutable std::vector<int>					mCurrentPositions;			///< Positions for the last codon decoded
};

#endif

