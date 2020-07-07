Aux étudiants en Master 2 de Bioinformatique: 
Bonjour, il est possible que vous, étudiant ayant pour projet de recoder le site protein interactions calculator, tombiez sur ce dépôt git (les sujet sont redondants d'une année à l'autre). Sachez que ce projet a été réalisé dans ce même objectif et dans le même contexte. J'ai décidé de laissé ce dépôt publique car il est un bon exemple de ce qu'il m'a été possible de concevoir en quelques jours seulement. Connaissant votre situation, je me doute que vous puissiez avoir envie de vous... inspirer de ce projet ;) . De ce fait, j'aimerais vous citer un petit passage de l'article L122-4 du code de la propriété intellectuelle : "toute représentation ou reproduction intégrale ou parielle faite sans le consentement de l'auteur ou de ses ayants droit ou ayants cause est illicite. Il en est de même pour la traduction, l'adaptation ou la transformation, l'arrangement ou la reproduction par un art ou un procédé quelconque."
Voila, à vous de voir si cela vaut le coup, en revanche, rien ne vous empeche de télécharger le script afin de comparer vos résultats (tips : il est possible / "normal" que vous ayez des sorties différentes du site, ayez un esprit critique).
Bon courage pour la réalisation de votre projet, et pour le reste du master, vous allez en avoir besoin !

P.S : venez au jebif, c'est sympa !

Bertrand Bouvarel, ancien M2BI, passionné d'analyse de données et de R


========== CREATION OF AN INTERACTION CALCULATION PROGRAM ==========
Author : Bertrand Bouvarel
Promotion : M2BI 2019-2020
DOI: 10.5281/zenodo.3442202 
last update: september 2019
Report available in the doc-report folder
Based on the article : Tina, K. G., Bhadra, R., & Srinivasan, N. (2007). 
		       PIC: Protein Interactions Calculator. Nucleic
                       Acids Research, 35(Web Server), W473–W476.

Python version : Python 3.7
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
intra and inter protein interactions, based on the atoms coordinates
of the PDB file of the protein.


========== USE ==========
Before running the program, please check that you have the main script
(cip.py) in in the src folder, as well as the packages (atom.py,
commandline.py, intcalc.py), in the packages folder. A file in PDB format
is required for launching the program. It is possible to access the help
of the program with the command (in the src folder):
>>>>> python3 cip.py --help <<<<<

In the src folder, run the program as follows in the bash:
>>>>> python3 cip.py ../data/file.pdb --intra --arg1 --arg2 --argN <<<<<

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
>>>>> python3 cip.py ../data/2c35.pdb --intra --inic --disu --arar --arsu --capi --hphb --mmhb --mshb --sshb <<<<<

>>>>> python3 cip.py ../data/2c35.pdb --inter --inic --disu --arar --arsu --capi --hphb --mmhb --mshb --sshb <<<<<

>>>>> python3 cip.py ../data/2c35.pdb --intra --inic6.5 --disu --arar4/6.5 --arsu5.5 --capi7 --hphb4.5 --mmhb --mshb --sshb <<<<<

>>>>> python3 cip.py ../data/2eti.pdb --intra --disu <<<<<
