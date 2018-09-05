
#ifndef NEWICKGRAMMAR_H
#define NEWICKGRAMMAR_H

#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_parse_tree.hpp>

using namespace boost::spirit::classic;

/// Grammar for the Newick formatted tree files.
/// Adapted from one of the Spirit examples
///
struct NewickGrammar : public grammar<NewickGrammar>
{
    static const int treeID      = 10;	///< Tree
    static const int nodelistID  = 20;	///< Nodelist
    static const int subtreeID   = 30;	///< Subtree
    static const int fulllabelID = 40;	///< Full label (name plus optional #1)
    static const int branchlenID = 50;	///< Branch length
    static const int cblenID     = 60;	///< Colon plus length
    static const int labelID     = 70;	///< Label
    static const int markerID    = 80;	///< Marker (i.e. the optional #1)

    template <typename ScannerT>
    struct definition
    {
        definition(NewickGrammar const& /*self*/)
        {
            tree = nodelist >> !full_label >> !colon_plus_len >> no_node_d[*space_p >> ch_p(';')];

            nodelist = no_node_d[ch_p('(')] >> subtree % (no_node_d[ch_p(',') >> *space_p]) >> no_node_d[ch_p(')')];

            subtree = nodelist >> !full_label >> !colon_plus_len >> no_node_d[*space_p] | full_label >> !colon_plus_len >> no_node_d[*space_p];

            colon_plus_len = no_node_d[ch_p(':') >> *space_p] >> branch_length;

            branch_length = real_p;

            full_label = no_node_d[ch_p('#')] >> marker | label >> !(no_node_d[ch_p('#')] >> marker);

            //label = leaf_node_d[alpha_p >> *(alnum_p|'.'|'_')];
            label = leaf_node_d[alnum_p >> *(alnum_p|'.'|'_')];

            marker = leaf_node_d[+alnum_p];
        }

        rule<ScannerT, parser_context<>, parser_tag<treeID> >      tree;
        rule<ScannerT, parser_context<>, parser_tag<nodelistID> >  nodelist;
        rule<ScannerT, parser_context<>, parser_tag<subtreeID> >   subtree;
        rule<ScannerT, parser_context<>, parser_tag<fulllabelID> > full_label;
        rule<ScannerT, parser_context<>, parser_tag<branchlenID> > branch_length;
        rule<ScannerT, parser_context<>, parser_tag<cblenID> >     colon_plus_len;
        rule<ScannerT, parser_context<>, parser_tag<labelID> >     label;
        rule<ScannerT, parser_context<>, parser_tag<markerID> >    marker;

        rule<ScannerT, parser_context<>, parser_tag<treeID> >   const&
        start() const { return tree; }
    };
};

#endif


