#!/bin/bash
echo "Digite o diretório de instalação:"
read DIR
cd ${DIR}
mkdir programs
cd programs
#
#Modulos Perl
apt-get install libgd-graph-perl
apt-get install libgraphviz-perl
apt-get install libxml-perl
apt-get install libxml-libxml-perl

perl -MCPAN -e 'install "Bio::SeqIO"'
perl -MCPAN -e 'install "Bio::AlignIO"'
perl -MCPAN -e 'install "Cwd"'
perl -MCPAN -e 'install "File::chdir"'
perl -MCPAN -e 'install "File::Copy"'
perl -MCPAN -e 'install "POSIX"'
perl -MCPAN -e 'install "Statistics::Distributions"'
perl -MCPAN -e 'install "Statistics::Multtest"'
perl -MCPAN -e 'install "Tie::File"'
perl -MCPAN -e 'install "Try::Tiny"'
perl -MCPAN -e 'install "Data::Dumper"'
perl -MCPAN -e 'install "File::Spec::Functions"'
perl -MCPAN -e 'install "File::Basename"'
perl -MCPAN -e 'install "FindBin"'
perl -MCPAN -e 'install "Capture::Tiny"'
perl -MCPAN -e 'install "Getopt::Long"'

#MUSCLE
mkdir muscle
cd muscle
wget -c "http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz"
tar -xzvf muscle3.8.31_i86linux64.tar.gz
cd ${DIR}/programs
#
#prank
mkdir prank
cd prank
wget -c "http://wasabiapp.org/download/prank/prank.linux64.140603.tgz"
tar -xzvf prank.linux64.140603.tgz
#echo "export PATH=$PATH:/opt/POTION-master/dp/prank/bin" >> ~/.bash_profile
echo "export PATH=$PATH:/opt/POTION_dev/POTION2/POTION2/programs/prank/prank/bin/* /usr/bin/"
cd ${DIR}/programs
#cp -R /opt/POTION-master/programs/prank/prank/bin/* usr/bin
#
#mafft
mkdir mafft
cd mafft
wget -c "https://mafft.cbrc.jp/alignment/software/mafft-7.313-gcc_fc6.x86_64.rpm"
apt-get install rpm
rpm -Uvh mafft-7.313-gcc_fc6.x86_64.rpm
cd ${DIR}/programs
#
#PhiPack
mkdir phipack
cd phipack
wget -c "www.maths.otago.ac.nz/~dbryant/software/PhiPack.tar"
tar -xf PhiPack.tar
cd PhiPack/src
make
cd ${DIR}/programs
#
#Phylip
mkdir phylip
cd phylip
wget -c "http://evolution.gs.washington.edu/phylip/download/phylip-3.697.tar.gz"
tar -xzvf phylip-3.697.tar.gz
cd phylip-3.697/src
cp Makefile.unx Makefile
make all
echo "export PATH=$PATH:`pwd`" >> ~/.bash_profile
cd ${DIR}/programs
#PALM
mkdir palm
cd palm
wget -c "http://abacus.gene.ucl.ac.uk/software/paml4.8a.macosx.tgz"
tar -xzvf paml4.8a.macosx.tgz
cd paml4.8
rm -f bin/*.exe
cd src
make -f Makefile
#mv baseml basemlg codeml pamp evolver yn00 chi2 ../bin
echo "export PATH=$PATH:${DIR}/programs/palm/paml4.8/bin/* /usr/bin/"
cd ${DIR}/programs
#TrimAl
mkdir trimal
cd trimal
wget -c "http://trimal.cgenomics.org/_media/trimal.v1.2rev59.tar.gz"
tar -xzvf trimal.v1.2rev59.tar.gz
cd trimAl/source
make
echo "export PATH=$PATH:/opt/POTION_dev/POTION2/POTION2/programs/trimal/trimAl/source/readal /usr/bin"
echo "export PATH=$PATH:/opt/POTION_dev/POTION2/POTION2/programs/trimal/trimAl/source/trimal /usr/bin"
cd ${DIR}/programs
#pagan
mkdir pagan
cd pagan
wget -c "http://wasabiapp.org/download/pagan/pagan.linux64.20150723.tgz"
tar -xzvf pagan.linux64.20150723.tgz
cd ${DIR}/programs
#phyml
mkdir phyml
cd phyml
wget -c "https://github.com/stephaneguindon/phyml/archive/master.zip"
unzip master.zip
cd phyml-master
./configure
make
cd ${DIR}/programs
#clustalo
mkdir clustalo
cd clustalo
apt-get install libargtable2-dev
wget -c "http://www.clustal.org/omega/clustal-omega-1.2.4.tar.gz"
tar -xzvf clustal-omega-1.2.4.tar.gz
cd clustal-omega-1.2.4/
./configure && make && make install
cd ${DIR}/programs
#Newick
mkdir newick
cd newick
wget -c "http://cegg.unige.ch/pub/newick-utils-1.6.tar.gz"
tar -xzvf newick-utils-1.6.tar.gz
cd newick-utils-1.6
./configure
make && make install



