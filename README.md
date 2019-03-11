# DeNovoCurvedSheetDesign
This repository contains the files necessary to generate most designs in the 2017 paper "Principles for designing proteins with cavities formed by curved beta sheets" (1)

<strong>Overview:</strong>

The code used to generate curved sheets following the principles is fully implemented in rosetta_scripts XMLs, with assisting python code to generate key input files: blueprints describing sheet register shifts and ABEGO bins, and constraint files dictating strand-strand hydrogen-bond pairing, curvature and twist. We have made working code for the folds described in the paper publicly available at this Git repository. The BPB_functions and Blueprint libraries contain a number of small functions used in the construction of all folds, while the code specific to each fold, used to prepare all the rosetta_script files, is located in each of the Fold_* folders, named PrepareFilesFold*.py. Each fold folder also comes with its own XML protocol, with the small variations necessary for each case.

<strong>Obtaining the Rosetta rosetta_scripts app:</strong>

The code in this repository generates all the files necessary for constructing the backbones using the application rosetta_scripts from the RosettaDesign software suit. 

Visit https://www.rosettacommons.org/docs/latest/Home for instructions on how to obtain and compile the latest version of rosetta_scripts.

<strong>Use example:</strong>

<code>$ cd Fold_D</code>

<code>$ ./PrepareFoldD.sh</code>

<code>$ ./run.sh </your/rosetta/exec/folder/rosetta_scripts.os.release> <runtime_in_sec></code>

<strong>General details of PrepareFilesFold*.py and Fold*_SheetGenerationProtocolTemplate.xml	files:</strong>

In each of the folders named after the folds designed in (1), there is a PrepareFilesFold*.py and a Fold*_SheetGenerationProtocolTemplate.xml file. PrepareFilesFold*.py generates blueprint and constraint files that direct the assembly process encoded in the Fold*_SheetGenerationProtocolTemplate.xml file. The backbone construction process is divided in at least four stages: 1) Construction of the two central strands of the main sheet. 2) Construction of the N-terminal flanking bulged strand. 3) Construction of the C-terminal flanking bulged strand. 4) Construction of the N-terminal helix/helices that pack agains the sheet. Step 4 is further divided in smaller steps for construction of more complex folds.

Each of the backbone construction steps is encoded by the same basic set of Rosetta movers and filters (see Fold*_SheetGenerationProtocolTemplate.xml), but uses different input files and values. This basic set is structured as a single mover that loops over constrain-guided fragment assembly followed by requirements checking, until all requirements are met or a number of attempts is reached. The constrain-guided fragment assembly is done in "centroid mode", which is basically a representation where side chains are approximated to a single soft sphere, with valine being the closest to the average size. This step also contains additional minimization steps to "pull" the structure together. The requirements checking step (or filtering step) contains a number of filters that change based on the step, but checks at least for Ramachandran outliers or unfavorable peptide bond geometry (e.g., cis peptide bonds). For specific documentation on movers and filtes, see the rosetta_scripts <a href="https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/RosettaScripts">wiki</a>.



<strong>Additional files:</strong>

Blueprint.py: Library by Javier Castellanos for Blueprint object handling.

BPB_functions.py: Library of functions for working with curved sheet protein backbones.

Rama_XPG_3level.txt: A "general" map from phi/psi coordinate to energy -  more populated Ramachandran space is more favorable. Author: Hahnbeom Park. All amino acids, except for proline and glycine, have the same values, hence the "general" qualifier. 

<strong>References:</strong>

(1) Marcos*, E., Basanta*, B., Chidyausiku, T. M., Tang, Y., Oberdorfer, G., Liu, G., ... & Pereira, J. H. (2017). Principles for designing proteins with cavities formed by curved β sheets. Science, 355(6321), 201-206.

*: Indicates co-first authorship
