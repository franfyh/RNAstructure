============
RNAstructure
============

This is my own version of RNAstructure, which includes tools for
RNA secondary structure prediction and analysis.

It is provided to you under the GNU GPL, which is spelled out
in gpl.txt.

In order to make the programs in the package run properly, please
change the ``DATAPATH`` environment variable.

In BASH, this is accomplished with:

``export DATAPATH=[directory in which RNAstructure resides]/RNAstructure/data_tables/``

In CSH, this is accomplished with:

``setenv DATAPATH [directory in which RNAstructure resides]/RNAstructure/data_tables/``

I contributed several parts to RNAstructure, they are:

Dynalign II
-----------

Under ``RNAstructure/`` directory, please run:

``make dynalign_ii``

The executable ``dynalign_ii`` will be in the ``RNAstructure/exe/`` directory.

The example configuration file to run Dynalign II on is
``RNAstructure/examples/dynalign_ii_sample.conf``. Please run:

``../exe/dynalign_ii dynalign_ii_sample.conf``

in ``RNAstructure/examples/`` directory.

Detailed instructions can be found at <http://rna.urmc.rochester.edu/Text/dynalign.html>.
Reference article can be found at <http://nar.oxfordjournals.org/content/42/22/13939>.

Multifind
---------
Under ``RNAstructure/`` directory, please run:

``make Multifind``

The executable ``Multifind`` will be in the ``RNAstructure/exe/`` directory.

The example configuration file to run Multifind on is 
``RNAstructure/examples/Multifind_sample.conf``

Please run:

``../exe/Multifind Multifind_sample.conf`` 

in ``RNAstructure/examples/`` directory.

Detailed instructions can be found at <http://rna.urmc.rochester.edu/Text/Multifind.html>.
Reference article can be found at <http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0130200>.

conserved_training
------------------
Under ``RNAstructure/`` directory, please run:

``make conserved_training``

The executable ``conserved_training`` will be in the ``RNAstructure/exe/`` directory.

The example configuration file to run conserved_training on is
``RNAstructure/conserved_training/conserved_training.conf``

Please run:

``../exe/conserved_training conserved_training.conf`` 

in ``RNAstructure/conserved_training/`` directory.

The ``TrainingFiles`` field of the conf file should include the names of the 
pairwise alignment
fasta file and the structure ct files for every data point in the training set.
The ``TestingFiles`` field of the conf file should include the names of the 
pairwise alignment
fasta file and the structure ct files for every data point in the testing set.
The ``Alpha`` field of the conf file is the parameter for the L2 regularization.

Reference article can be found at <https://github.com/franfyh/RNAstructure/blob/master/maximum_entropy.pdf>.

TurboFold
---------
The improvement of TurboFold is a project ongoing. Original TurboFold can be 
downloaded at <http://rna.urmc.rochester.edu/RNAstructure.html>.
The reference of the original TurboFold can be found at <http://www.biomedcentral.com/1471-2105/12/108>.
