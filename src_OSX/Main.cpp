// **************************************************************************************************
// biLouvain: A C++ library for bipartite graph community detection
// Paola Gabriela Pesantez-Cabrera, Ananth Kalyanaraman  
//	(p.pesantezcabrera@wsu.edu, ananth@eecs.wsu.edu)
// Washington State University
//
// For citation, please cite the following paper:
// Pesantez, Paola and Kalyanaraman, Ananth, "Detecting Communities in Biological 
// Bipartite Networks," Proc. ACM Conference on Bioinformatics, Computational Biology, 
// and Health Informatics (ACM-BCB), Seattle, WA, October 2-5, 2016, In press, 
// DOI: http://dx.doi.org/10.1145/2975167.2975177.
//
// **************************************************************************************************
// Copyright (c) 2016. Washington State University ("WSU"). All Rights Reserved. 
// Permission to use, copy, modify, and distribute this software and its documentation 
// for educational, research, and not-for-profit purposes, without fee, is hereby 
// granted, provided that the above copyright notice, this paragraph and the following 
// two paragraphs appear in all copies, modifications, and distributions. For 
// commercial licensing opportunities, please contact The Office of Commercialization, 
// WSU, 280/286 Lighty, PB Box 641060, Pullman, WA 99164, (509) 335-5526, 
// commercialization@wsu.edu<mailto:commercialization@wsu.edu>, https://commercialization.wsu.edu/

// IN NO EVENT SHALL WSU BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, 
// OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF 
// THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF WSU HAS BEEN ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE.

// WSU SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND 
// ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". WSU HAS NO 
// OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
// **************************************************************************************************

#include "Header.h"
#include "PreProcessInputBipartiteGraph.h"
#include "Graph.h"
#include "LoadGraph.h"
#include "LoadSimilarityMatrix.h"
#include "FuseMethod.h"
#include "biLouvainMethod.h"
#include "biLouvainMethodMurataPN.h"
#include "Timer.h"


static std::string inputFileName = "";
static std::string initialCommunitiesFileName = "";
static std::string similarityMatrixFileName = "";
static std::string outputFileName = "";
static std::string delimiter = "\t";
static int optionOrder = 3;
static int fuse = 1;
static double cutoffIterations = 0.01;
static double cutoffPhases =  0.0;
static double cutoffFuse = 1.0;
static double alpha =  1.0;
static int flag;
static void parseCommandLine(const int argc, char * const argv[]);

static struct option longopts[] = {
   { "input",		required_argument,0,'i'},
   { "delimeter",	required_argument,0,'d'},
   { "ci",		required_argument,&flag,1},
   { "cp",		required_argument,&flag,2},
   { "order",		required_argument,&flag,3},
   { "initial",		required_argument,&flag,4},
   { "fuse",		required_argument,&flag,5},
   { "cf",              required_argument,&flag,6},
   { "similarity",	required_argument,&flag,7},
   { "alpha",		required_argument,&flag,8},
   { "output",		required_argument,0,'o'},
   { 0, 0, 0, 0 }
};


int main(int argc, char *argv[])
{
	try
	{
		parseCommandLine(argc, argv);
		struct timeval startTime,endTime;	
		std::ifstream infile(inputFileName.c_str());
		//std::cout<<inputFileName<<std::endl;
		int items = 0;
                bool band = false;
                if(infile.is_open()==true)
                {
                        std::cout << "\n ::: Loading Bipartite Graph " << inputFileName << " :::";
                        std::string bipartiteFileName;
                        std::unordered_map<int,std::string> bipartiteOriginalEntities;
                        //if((inputFileName.find("bipartite") == std::string::npos)||(inputFileName.find("Bipartite") == std::string::npos))
                        std::string* pieces;
                        pieces = StringSplitter::split(inputFileName,"bipartite",items);
                        if(items > 1)
                                band = true;
                        else
                        {
                                pieces = StringSplitter::split(inputFileName,"Bipartite",items);
                                if(items > 1)
                                        band = true;
				else
				{
					pieces = StringSplitter::split(inputFileName,"BIPARTITE",items);
                                	if(items > 1)
                                        	band = true;
				}
                        }
                        delete[] pieces;
                        if(band == false)
                        {
				int pos = inputFileName.find_last_of(".");
				bipartiteFileName = inputFileName.substr(0,pos)+"_bipartite.txt";
				bipartiteOriginalEntities=PreProcessInputBipartiteGraph::preProcessingGraphData(inputFileName,delimiter);
			}
			else
			{
				bipartiteFileName = inputFileName;
				bipartiteOriginalEntities=PreProcessInputBipartiteGraph::readDictionaryFile(inputFileName);
			}
			double loadGraphTime = 0.0;
			int pass = -1;
			gettimeofday(&startTime,NULL);
			Graph* graph;
			pass = LoadGraph::loadBipartiteGraphFromFile(graph,bipartiteFileName);
			gettimeofday(&endTime,NULL);
			loadGraphTime = (endTime.tv_sec - startTime.tv_sec)*1000000 + (endTime.tv_usec - startTime.tv_usec);
			if (pass == 0)
			{                            
				std::cout << "\n ::: Done Loading Bipartite Graph :::";		
				if(alpha != 1.0)
				{
					int numberMatrices = 0;
					std::string similarityMatrixPartitionFileName = "";
					if(similarityMatrixFileName.empty()== false)
					{	
						similarityMatrixPartitionFileName = similarityMatrixFileName + "V1.csv";
						pass = LoadSimilarityMatrix::loadSimilarityMatrixFromFile(*graph,similarityMatrixPartitionFileName); 
						if (pass == 0)	
							numberMatrices ++;
						else
							printf("\n ::: Warning: Similarity Matrix for vertices in V1 was not found. :::");
						similarityMatrixPartitionFileName = similarityMatrixFileName + "V2.csv";
						pass = LoadSimilarityMatrix::loadSimilarityMatrixFromFile(*graph,similarityMatrixPartitionFileName); 
						if (pass == 0)	
							numberMatrices ++;
						else
							printf("\n ::: Warning: Similarity Matrix for vertices in V2 was not found. :::");
					}	
					//std::cout << "\n ::: Alpha V1: " << graph->getAlphaV1() << "\t Alpha V2: " << graph->getAlphaV2();
					if(numberMatrices == 0)
						alpha = 1.0;
				}	
				FuseMethod f;
				biLouvainMethodMurataPN biLouvain;
				//std::vector<double> communitiesBetaFactor;
				 if((fuse == 1)&&(initialCommunitiesFileName.empty()==true))
                                        f.fuseMethodFile(*graph,bipartiteFileName,alpha,cutoffFuse);
                                else if((fuse == 1)&&(initialCommunitiesFileName.empty()==false))
                                        f.initialCommunityDefinitionProvidedFileCommunities(*graph,initialCommunitiesFileName,alpha);
				std::cout << "\n ::: Starting biLouvain Algorithm :::";
				gettimeofday(&startTime,NULL);							
				if(alpha != 0.0)
					biLouvain.biLouvainMethodAlgorithm(*graph,cutoffIterations,cutoffPhases,optionOrder,bipartiteOriginalEntities,bipartiteFileName,outputFileName,alpha);
				else
					biLouvain.biLouvainMethodAlgorithmIntraType(*graph,cutoffIterations,cutoffPhases,optionOrder,bipartiteOriginalEntities,bipartiteFileName,outputFileName);
				gettimeofday(&endTime,NULL);	
				double biLouvainAlgorithmTime = (endTime.tv_sec - startTime.tv_sec)*1000000 + (endTime.tv_usec - startTime.tv_usec);
				biLouvain.printTimes(biLouvainAlgorithmTime,loadGraphTime,f.fusingTime);
				graph->destroyGraph();
			}
			else
			{
				printf("\n ::: There was a problem reading the graph input file :::\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
				printf("\n ::: Input file was not found :::\n");
				exit(EXIT_FAILURE);
		}
	}
	catch (...)
	{
		std::cout << " ::: Unknown error Main :::" << std::endl;
	}
	std::cout << std::endl << " ::: biLouvain Method has finished :::" << std::endl;
	return 0;
}

void printUsage()
{
	 printf("Usage: -i {inputFile} -d {delimiter (\" \",\",\",\"\\t\")} [-ci {cutoff iterations(default=0.01)} -cp {cutoff phases(default=0.0)} -order {1:Sequential, 2:Alternate, 3:Random(default=3)} -initial {initialCommuitiesFile}(default=\"\") -fuse {0/1 flag(default=1)} -o {outputFileName(default=input_Results*)}]\n");  
         exit(EXIT_FAILURE);
}


void parseCommandLine(const int argc, char * const argv[])
{
	int c=0,indexPtr=0,prevInd;
	while(prevInd=optind,(c = getopt_long_only(argc, argv, "i:d:o:", longopts, &indexPtr)) != -1) {
		//printf("%d \t %d \t %c\n",prevInd,optind,c);
		if((optind==prevInd+2) && (*optarg=='-'))
		{
			printf("::: You forgot to provide an argument %s :::\n");
			printUsage();
		}
		switch (c) {
		    case 'i':
		        inputFileName = optarg;
			if(inputFileName.empty())
			{	printf(" ::: No input file name provided :::\n");
				printUsage();
			}
		    	break;
		    case 'd':
		        delimiter = optarg;
			if(delimiter.empty())
			{	printf(" ::: No delimiter provided :::\n");	
				printUsage();
			}
			else if((delimiter != " ")&&(delimiter != "\\t")&&(delimiter != ","))
			{
				printf(" ::: Unknown delimiter provided :::\n");
		        printUsage();		
			}
			else if(delimiter == "\\t")
				delimiter = "\t";
			break;
		    case 'o':
			if(optarg != NULL)
				outputFileName = optarg;
			break;	
		    case 0:
			if(*(longopts[indexPtr].flag)==1)
			{
				if(optarg != NULL)
		                      cutoffIterations = atof(optarg);
			}
			else if(*(longopts[indexPtr].flag)==2)
			{
				if(optarg != NULL)
					cutoffPhases = atof(optarg);
			}
			else if(*(longopts[indexPtr].flag)==3)
			{
				if(optarg != NULL)
					optionOrder = atoi(optarg);
			}
			else if(*(longopts[indexPtr].flag)==4)
			{
				if(optarg != NULL)
					initialCommunitiesFileName = optarg;
			}
			else if(*(longopts[indexPtr].flag)==5)
			{
				if(optarg != NULL)
					fuse = atoi(optarg);
			}
			else if(*(longopts[indexPtr].flag)==6)
                        {
                                if(optarg != NULL)
                                        cutoffFuse = atof(optarg);
                        }
			else if(*(longopts[indexPtr].flag)==7)
			{
				if(optarg != NULL)
					similarityMatrixFileName = optarg;
			}
			else if(*(longopts[indexPtr].flag)==8)
			{
				if(optarg != NULL)
					alpha = atof(optarg);
			}
			break;
		    case ':':
			printUsage;
			break;
		    case '?':
			/*if((optopt == 'i') || (optopt == 'd'))
				printf("Argument is mandatory for -%c\n",optopt);
			else if(isprint(optopt))
				printf("You have provided an unknown option -%c\n",optopt);
			else
				printf("Unknown Error-0x%08x\n",optopt);			
			*/
			printUsage();
			break;
		    default:   
			exit(0);
		}
	}
}
