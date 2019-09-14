========== CREATION OF AN INTERACTION CALCULATION PROGRAM ==========
Author : Bertrand Bouvarel
Promotion : M2BI 2019-2020
Python version : Python 3.7.4
Based on the article : Tina, K. G., Bhadra, R., & Srinivasan, N. (2007). 
		       PIC: Protein Interactions Calculator. Nucleic
                       Acids Research, 35(Web Server), W473–W476.

Packages used :
- sys
- math
- statistics
- command_line (written by the author)
- atom (written by the author)
- intcalc (written by the author)


========== DESCRIPTION ==========
Program for the calculation of disulphide bridges, hydrogen bonds and
hydrophobic, ionic, aromatic-aromatic, aromatic-sulfur, cation-pi
interactions, based on the atoms coordinates of the PDB file of the
protein.


========== USE ==========
Before running the program, please check that you have the main script
(cip.py) in in the bin folder, as well as the packages (atom.py,
commandline.py, intcalc.py), in the packages folder. A file in PDB format
is required for launching the program. It is possible to access the help
of the program with the command (in the bin folder):
>>>>> python3 cip.py --help <<<<<

In the bin folder, run the program as follows in the bash:
>>>>> python3 cip.py ../data/file.pdb --arg1 --arg2 --argN <<<<<

Will be displayed at the end of the execution and written in a file 
named as the protein, in the results folder, the list of the interactions
of the protein, calculated with the coordinates of the atoms of the
protein. Will be displayed in these tables:
- the two positions of the residues
- the residues doing the interaction
- the chain of the residues
- the name of the atom involved in the interaction (hydrogen bonds only)
- the distance between the atoms / centroids involved in the interaction


========== EXAMPLES ==========
To launch the calculation of all of the interactions implemented with 
the default values of distance, use the command:
>>>>> python3 cip.py ../data/1atn.pdb --inic --disu --arar --arsu --capi --hphb --mmhb --mshb --sshb <<<<<

>>>>> python3 cip.py ../data/1bta.pdb --inic --disu --arar --arsu --capi --hphb --mmhb --mshb --sshb <<<<<

>>>>> python3 cip.py ../data/2eti.pdb --inic --disu --arar --arsu --capi --hphb --mmhb --mshb --sshb <<<<<

