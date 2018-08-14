# DeNovoCurvedSheetDesign
This repository contains the files necessary to generate most designs in the 2017 paper "Principles for designing proteins with cavities formed by curved beta sheets"

The code used to generate curved sheets following the principles is fully implemented in rosetta_scripts XMLs, with assisting python code to generate key input files: blueprints describing sheet register shifts and ABEGO bins, and constraint files dictating strand-strand hydrogen-bond pairing, curvature and twist. We have made working code for the folds described in the paper publically available at this Git repository. The BPB_functions and Blueprint libraries contain a number of small functions used in the construction of all folds, while the code specific to each fold, used to prepare all the rosetta_script files, is located in each of the Fold_* folders, named PrepareFilesFold*.py. Each fold folder also comes with its own XML protocol, with the small variations necessary for each case.

Blueprint.py: Library by Javier Castellanos for Blueprint object handling.

BPB_functions.py: Library of functions for working with curved sheet protein backbones.

Rama_XPG_3level.txt: A "general" map from phi/psi coordinate to energy -  more populated Ramachandran space is more favorable. All amino acids, except for proline and glycine, have the same values, hence the "general" qualifier. 
