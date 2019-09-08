==========CREATION D'UN PROGRAMME DE FLIP D'ASPARAGINES ET DE GLUTAMINES==========
Projet proposé par : Romain Retureau
Contributeurs : Bertrand Bouvarel, Alexandre Lecoeur
Promotion : M1BI 2018-2019
Version de Python : Python 3.7

Modules utilisés :
- math
- sys


==========DESCRIPTION SIMPLE==========
Ce programme permet, à partir d'un fichier pdb et d'un fichier topologique de sortir 
une nouveau fichier pdb avec les modification des position des asparagines et glutamines 
suite aux flip permetant d'obtenir une énergie électrostatique plus basse.


==========UTILISATION==========	
Avant de lancer le programme, merci de vérifier que vous possédez bien le script 
ainsi que le fichier pdb de la séquence (nucleosome_3MVD_1KX5_noH.pdb)et le fichier 
topologique (toppar/top_all36_prot.rtf) dans le même répertoire. Le fichier topologique
restera le même peu importe la séquence donc le fichier pdb utilisé.
Dans le répertoire contenant ces éléments, lancer le programme de la manière suivante
dans le bash:
>>>>>>>>>> python3 script.py nucleosome_3MVD_1KX5_noH.pdb top_all36_prot.rtf <<<<<<<<<<

S'affichera à la fin de l'éxecution:
Le nombre de flip pour chacun des deux acides aminés afin de se rendre compte de 
l'ensemble des modifications réalisées par rapport à la séquence initiale.
