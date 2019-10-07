#!/bin/bash -l
#SBATCH -J biLouvainMethod
#SBATCH --mail-user=p.pesantezcabrera@wsu.edu
#SBATCH -q debug
#SBATCH -N 1
#SBATCH -t 00:10:00
#SBATCH --mem 10GB
#SBATCH -L project
#SBATCH --mail-type=ALL


module load gcc/6.1.0
#folder="/global/homes/p/ppesante/biLouvain/inputData/ThesisResults/fuse0bin/"
folder="/global/homes/p/ppesante/biLouvain/inputData/"
#folder="/global/homes/p/ppesante/biLouvain/inputData/interactions/Base/"
executable=biLouvain
numProcess=(1)
numThreads=(1)
order=3
cutoffPhases=0.0
cutoffIterations=0.01

#datasets=(AuthorPaper scotland)
#datasets=(Safariland barrett1987 bezerra2009 elberling1999 inouye1988 junker2013 kato1990 kevan1970 memmott1999 mosquin1967 motten1982 olesen2002aigrettes olesen2002flores ollerton2003 schemske1978 small1976 vazarr vazcer vazllao vazmasc vazmasnc vazquec vazquenc)
#datasets=(SouthernWomen testCaseID DivorceUS scotland host-rolesappID AuthorPaper malaria complexesDrug gene-drug-ID)
#datasets=(SouthernWomen testCaseID DivorceUS scotland host-rolesappID AuthorPaper malaria complexesDrug gene-drug-ID UpE115L3ID UpE135L3ID UpE155L3ID UpP4L3ID UpP28L3ID UpP14L3ID UpE185L3ID edgeListBig eQTL ebola)
#datasets=(GeneDisease2016-ID)
#datasets=(SpeciesInteractions_EID2ID) 
#datasets=(UpE115L3ID UpE135L3ID UpE155L3ID UpP4L3ID UpP28L3ID UpP14L3ID UpE185L3ID) 
#datasets=(UpE115L3ID UpE155L3ID UpP4L3ID UpP28L3ID)
#datasets=(UpP14L3ID UpP28L3ID SpeciesInteractions_EID2ID)
#datasets=(junker2013 kato1990 kevan1970 memmott1999 malaria complexesDrug)
#datasets=(edgeList_wikiPW)
#datasets=(enzime_interactions gpcr_interactions ion_chanel_interactions nr_interactions)
#datasets=(LPBrimBipartitekato1990 LPBrimBipartitekevan1970)
datasets=(biSBM301 biSBM30)
#datasets=(BCRing22 BCRing23 BCRing33 BCRing34 BCRing56 BCRingM22 BCRingM56 BCRingI22 BCRingI56 BCRing1020 BCChain203 BCChain264Com BCChain550 BCChain6100)
#datasets=(eQTL testCaseID DivorceUS scotland host-rolesappID)
for i in ${datasets[@];}
do
        fileName=$folder${i}"_bipartite.txt"
	initialCommunities=$folder${i}"_InitialCommunities.txt"
        srun -n $numProcess ./$executable -i $fileName -d "," -order $order -ci $cutoffIterations -cp $cutoffPhases -fuse 1
#       srun -n $numThreads valgrind --tool=memcheck --leak-check=yes --track-origins=yes ./$executable $fileName $whatToDo $order $cutoffIterations $cutoffPhases $dictionaryName
done
exit 0

