// **************************************************************************************************
// biLouvain: A C++ library for bipartite graph community detection
// Paola Gabriela Pesantez-Cabrera, Ananth Kalyanaraman
//      (p.pesantezcabrera@wsu.edu, ananth@eecs.wsu.edu)
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



#include "FuseMethod.h"

FuseMethod::FuseMethod():biLouvainMethodMurataPN()
{
	fusingTime = 0.0;
}
FuseMethod::~FuseMethod(){}

int FuseMethod::fuseCommunities(Graph &g,int start,int end, double &lambda)
{
	int maxIntersection = 0;
	int result = 0;
        int key = 0;
        std::vector<int> a;
        std::vector<int> b;
	std::vector<int> c;
	std::vector<int> coClusterMate;
	double temp = 0.0, betaFactor = 0.0;
	double totalModularity = 0.0;
	int first = 0;
	if(_alpha != 0.0)
		totalModularity = CoClusterMateDefinitionAllCommunities(g,start,end);
	else
		totalModularity = IntraTypeDefinitionAllCommunities(g,start,end);
	for(int i=start;i<end;i++)
	{
	  maxIntersection = 0;
	  betaFactor = 0.0;
          key = i;
	  if(_alpha != 0.0)
	         a =  _communities[i].getCoClusterMateCommunityId();
          //std::unordered_set<int> c(a.begin(),a.end());
          for(int j=start;j<i;j++)
          {
         	//std::cout << "\nCommunity: " << j << "\t Number Nodes: " << _communities[j].getNumberNodes() << std::endl;
                if(_communities[j].getNumberNodes()>0) 
		{
                	if(_alpha != 0.0)
				b = _communities[j].getCoClusterMateCommunityId();
                        //int intersection = std::count_if(b.begin(),b.end(),[&](int k){return c.find(k) != c.end();});
			//similarityCommunity = calculateCommunitySimilarity(i);	
			//std::cout << "\nBeta Factor entre: " << key << "\t"<<j<<"\t"<< temp;
			if((_alpha > 0.0)&&(_alpha < 1.0))
                        {
				set_intersection(a.begin(),a.end(),b.begin(),b.end(),back_inserter(c));
				temp = calculateCommunityBetaFactor(g,_communities[j].getDescription(),g._graph[i].getSimilarityToIntraTypeNeighborCommunity(j));
				if((c.size() > maxIntersection)&&(temp >= lambda)&&(temp>betaFactor))
                        	{
                        		key = j;
                                	maxIntersection = c.size();
					coClusterMate.clear();
					coClusterMate = c;
					betaFactor = temp;
				}
                        }
			else if(_alpha == 1.0)
			{
				set_intersection(a.begin(),a.end(),b.begin(),b.end(),back_inserter(c));
				if(c.size() > maxIntersection)
                                {
                                        key = j;
                                        maxIntersection = c.size();
                                        coClusterMate.clear();
                                        coClusterMate = c;
                                }
			}
			else
			{
				temp = calculateCommunityBetaFactor(g,_communities[j].getDescription(),g._graph[i].getSimilarityToIntraTypeNeighborCommunity(j));
				if((temp >= lambda)&&(temp > betaFactor))
                                {
                                        key = j;
                                        betaFactor = temp;
					//std::cout << "\nCommunity: " << j << "\tBF: " << betaFactor<< "\tTemp: "<< temp << std::endl;
                                }
			}
                        b.clear();
			c.clear();
                        //std::cout << "\nCommunity Comparison: " << j << "\t Inter: " << intersection << std::endl;
                }
         }
         a.clear();
         if(key != i) //we have communities to merge
         {
		result++;
                //update communities
		updateNodeCommunity(g,g._graph[i].getId(),i,key);
		updateNeighborCommunities(g,g._graph[i].getId(),i,key);
		if(_alpha != 0.0)
			_communities[i].setCoClusterMateCommunityId(coClusterMate);
		if(_alpha != 1.0)
                {
			updateNodeIntraTypeCommunity(g,g._graph[i].getId(),i,key);
                        updateIntraTypeNeighborCommunities(g,g._graph[i].getId(),i,key);
			_communities[i].setBetaFactor(0.0);
	                _communities[key].setBetaFactor(betaFactor);	
			//printf("\nPrint del Fuse::: Community: %d    BF: %f",key,betaFactor);
                }
         }
         //std::cout << "\nCommunity: " << i << "   Key:  " << key << std::endl;
	}
	return result;
}

double FuseMethod::fuseMethodInit(Graph &g)
{
	double finalModularity = 0.0;
	if(_alpha == 1.0)
        {
                initialCommunityDefinition(g);
                initialCommunityNeighborsDefinition(g);
        }
        else
        {
                initialCommunityDefinitionWithIntraType(g);
                initialCommunityNeighborsDefinition(g);
                initialIntraTypeCommunityNeighborsDefinition(g);
        }
        if(_alpha != 0.0)
                finalModularity = CoClusterMateDefinitionAllCommunities(g,0,_numberCommunities);
        else
                finalModularity = IntraTypeDefinitionAllCommunities(g,0,_numberCommunities);
	return finalModularity;
}

void FuseMethod::fuseMethodCalculationMF(Graph &g, std::string outputFileName,double cutoffFuse)
{
	double finalModularity = 0.0, initialModularity = 0.0;
	int changes = 0;
	struct timeval startTime,endTime;
        double fuseTime = 0.0;
        finalModularity = fuseMethodInit(g);
        do
        {
                gettimeofday(&startTime,NULL);
                initialModularity = finalModularity;
                changes = fuseCommunities(g,0,g._lastIdPartitionV1+1,g._lambdaV1);
                changes += fuseCommunities(g,g._lastIdPartitionV1+1,_numberCommunities,g._lambdaV2);
                if(changes > 0)
                {
                        fromCommunitiesToNodes(g);
			finalModularity = fuseMethodInit(g);
                        gettimeofday(&endTime,NULL);
                        fuseTime = (endTime.tv_sec - startTime.tv_sec)*1000000 + (endTime.tv_usec - startTime.tv_usec);
                        std::cout<<"\n ::: Fuse "<<finalModularity<<"\t"<<initialModularity<<"\t"<<g._lastIdPartitionV1+1<<"\t"<<g._numberNodes-(g._lastIdPartitionV1+1)<<"\t"<<g._numberNodes<<"\t"<<g._numberEdges<<"\t"<<fuseTime<<std::endl;
                }
        }while((finalModularity - initialModularity) > cutoffFuse);
	
	std::stringstream initialCommunities;
        std::ofstream outputFile;
        outputFile.open(outputFileName.c_str(),std::ios::out|std::ios::trunc);
        std::vector<int>community;
        for(int i=0;i<_numberCommunities;i++)
        {
        	if(_communities[i].getNumberNodes()>0)
                {
                	for(int j=0;j<g._graph[_communities[i].getNodes()[0]].getNumberNodes();j++)
                        	initialCommunities << g._graph[_communities[i].getNodes()[0]].getNodesSorted()[j].getIdInput() << ",";
                        initialCommunities.seekp(initialCommunities.str().length()-1);
                        initialCommunities << "\n";
                        //std::cout<<initialCommunities.str();
                        outputFile << initialCommunities.str();
                        initialCommunities.str("");
                        community.clear();
                }
         }
         outputFile.close();
	//for(int i=0;i<g._numberNodes;i++)
        //        std::cout << g._graph[i].getId()<< "\t" << g._graph[i].getDegreeNode()<<std::endl;
	//return communitiesBetaFactor;
}

void FuseMethod::fuseMethodCalculation(Graph &g, std::string outputFileName)
{
	int changes = 0;
	if(_alpha == 1.0)
	{
		initialCommunityDefinition(g);                                                 
              	initialCommunityNeighborsDefinition(g);
	}
	else
	{
		initialCommunityDefinitionWithIntraType(g);                                               
                initialCommunityNeighborsDefinition(g);
		initialIntraTypeCommunityNeighborsDefinition(g);
	}	
	//double totalModularity = CoClusterMateDefinitionAllCommunities(g,0,_numberCommunities,alpha); //assign the cocluster community(ies) to each community
	//std::cout << "\nInitial Modularity: " << totalModularity << std::endl; 
	changes = fuseCommunities(g,0,g._lastIdPartitionV1+1,g._lambdaV1);
	changes += fuseCommunities(g,g._lastIdPartitionV1+1,_numberCommunities,g._lambdaV2);
	std::stringstream initialCommunities;
        std::ofstream outputFile;
        outputFile.open(outputFileName.c_str(),std::ios::out|std::ios::trunc);
	std::vector<int>community;
	for(int i=0;i<_numberCommunities;i++)
	{
		if(_communities[i].getNumberNodes()>0)
		{
			community = _communities[i].getNodes();
			sort(community.begin(),community.end());
			for(int j=0;j<_communities[i].getNumberNodes();j++)
				initialCommunities << community[j] << ",";
	                initialCommunities.seekp(initialCommunities.str().length()-1);
			initialCommunities << "\n";
	                //std::cout<<initialCommunities.str();
        	        outputFile << initialCommunities.str();
                	initialCommunities.str("");
			community.clear();
		}
        }
	//printCommunitiesContributionModularity();
        outputFile.close();
	fromCommunitiesToNodes(g);
	std::cout<<"\n ::: Fuse "<<g._lastIdPartitionV1+1<<"\t"<<g._numberNodes-(g._lastIdPartitionV1+1)<<"\t"<<g._numberNodes<<std::endl;
	//for(int i=0;i<g._numberNodes;i++)
        //        std::cout << g._graph[i].getId()<< "\t" << g._graph[i].getDegreeNode()<<std::endl;
}

void FuseMethod::initialCommunityDefinitionProvidedFileCommunities(Graph &g,const std::string &initialCommunitiesFileName,double &alpha)
{
	_alpha = alpha;
        std::ifstream initialCommunitiesFile(initialCommunitiesFileName.c_str());
        std::unordered_map<int,double> nodesInCommunity;
	std::unordered_map<int,double> nodesIntraType;
	//std::vector<double> communitiesBetaFactor;
        int numberCommunities = 0;
        int numberNodes = 0;
        int id = 0;
        std::string line;
        std::string* nodes;
	std::string* lineCommunity;
	int items = 0;
        if(initialCommunitiesFile.is_open())
        {
		while(initialCommunitiesFile.good())
                {
                        getline(initialCommunitiesFile,line);
                        if(initialCommunitiesFile.eof())break;
                        if(line.find(",")!= std::string::npos)
                               	nodes = StringSplitter::split(line,",",numberNodes);
                        else
                               	nodes = StringSplitter::split(line,"\n",numberNodes);
			for(int i=0;i<numberNodes;i++)
                        {
				//std::cout<<nodes[i]<<",";
                                id = stoi(nodes[i]);
				g._graph[id].setCommunityId(numberCommunities);
	                        nodesInCommunity[g._graph[id].getId()] = g._graph[id].getDegreeNode(); 
				if(_alpha != 1.0)
					nodesIntraType[g._graph[id].getId()] = g._graph[id].getSimilarityNode();
			}
			//std::cout<<"\n";
                        Community community(numberCommunities,g._graph[id].getType(),nodesInCommunity,nodesIntraType);
                        _communities.push_back(community);
			/*for(int j=0;j<_communities[numberCommunities].getNumberNodes();j++)
                                std::cout << _communities[numberCommunities].getNodes()[j] << ",";
                        std::cout<<"\n";*/
                        nodesInCommunity.clear();
			nodesIntraType.clear();
                        numberCommunities++;
                }
                _numberCommunities = numberCommunities;
                delete[] nodes;
                initialCommunitiesFile.close();
		fromCommunitiesToNodes(g);
		std::cout<<"\n ::: Fuse "<<g._lastIdPartitionV1+1<<"\t"<<g._numberNodes-(g._lastIdPartitionV1+1)<<"\t"<<g._numberNodes<<std::endl;
	        //for(int i=0;i<g._numberNodes;i++)
                //      std::cout << g._graph[i].getId()<< "\t" << g._graph[i].getDegreeNode()<<std::endl;
		//return communitiesBetaFactor;
        }
        else
        {
                printf("\nInitial Communities File not found\n");
                exit(EXIT_FAILURE);
        }
}

void FuseMethod::initialCommunityDefinitionProvidedFileMetaNodes(Graph &g,const std::string &initialCommunitiesFileName,double &alpha)
{
	_alpha = alpha;
        std::ifstream initialCommunitiesFile(initialCommunitiesFileName.c_str());
	int lastIdPartitionV1 = -1;
        std::vector<Node> nodesInCommunity;
        std::unordered_map<int,double> neighbors;
	std::vector<MetaNode>newGraph;
	int numberCommunities = 0;
        int numberNodes = 0;
        int id = 0;
        std::string line;
        std::string* nodes;
        if(initialCommunitiesFile.is_open())
        {
                while(initialCommunitiesFile.good())
                {
                        getline(initialCommunitiesFile,line);
                        if(initialCommunitiesFile.eof())break;
                        if(line.find(",")!= std::string::npos)
                                nodes = StringSplitter::split(line,",",numberNodes);
                        else
                                nodes = StringSplitter::split(line,"\n",numberNodes);
			for(int i=0;i<numberNodes;i++)
			{	
				//std::cout<<nodes[i]<<",";
				id = stoi(nodes[i]);
				g._graph[id].setCommunityId(numberCommunities);
                                nodesInCommunity.push_back(g._graph[id].getNodes()[0]);
                        }
			//std::cout<<"\n";
			MetaNode metanode(numberCommunities,nodesInCommunity[0].getType(),nodesInCommunity,neighbors,-1);
                        numberCommunities++;
			newGraph.push_back(metanode);
                        if(metanode.getType()=="V1")
                                lastIdPartitionV1++;
                        nodesInCommunity.clear();
                }
		for(int i=0;i<numberCommunities;i++)
		{
			for(int j=0;j<newGraph[i].getNumberNodes();j++)
			{
				for(int k=0; k<g._graph[newGraph[i].getNodes()[j].getIdInput()].getNumberNeighbors();k++)
                        	{
                        		int idNeighbor = g._graph[g._graph[newGraph[i].getNodes()[j].getIdInput()].getNeighbors()[k]].getCommunityId();
	                                if(neighbors.find(idNeighbor)!= neighbors.end())
        	                        	neighbors[idNeighbor] += g._graph[newGraph[i].getNodes()[j].getIdInput()].getWeightNeighbor(g._graph[newGraph[i].getNodes()[j].getIdInput()].getNeighbors()[k]);
                	                else
                                	        neighbors[idNeighbor] = g._graph[newGraph[i].getNodes()[j].getIdInput()].getWeightNeighbor(g._graph[newGraph[i].getNodes()[j].getIdInput()].getNeighbors()[k]);
                        	        //std::cout << "Neighbor:" << idNeighbor << "  Weight: " << neighbors[idNeighbor] << std::endl;
				}
                        }
			newGraph[i].setNeighbors(neighbors);
			neighbors.clear();
		}
                delete[] nodes;
                initialCommunitiesFile.close();
		MetaNode*_newGraph = new MetaNode[numberCommunities];
		for(int i=0;i<numberCommunities;i++)
			_newGraph[i] = newGraph[i];
		Graph compactedGraph(_newGraph,numberCommunities,g._numberEdges,g._weightEdges,g._weightEdgesV1,g._weightEdgesV2,lastIdPartitionV1);
		g.destroyGraph();
	        g = compactedGraph;
		std::cout<<g._lastIdPartitionV1+1<<"\t"<<g._numberNodes-(g._lastIdPartitionV1+1)<<"\t"<<g._numberNodes<<std::endl;
                //for(int i=0;i<g._numberNodes;i++)
                //      std::cout << g._graph[i].getId()<< "\t" << g._graph[i].getDegreeNode()<<std::endl;
        }
        else
        {
                printf("\nInitial Communities File not found\n");
                exit(EXIT_FAILURE);
        }
}

void FuseMethod::fuseMethodFile(Graph &g,const std::string &inputFileName,double &alpha, double cf)
{
	//std::vector<double> communitiesBetaFactor;
	_alpha = alpha;
	struct timeval startTime,endTime;		
	int pos = inputFileName.find_last_of(".");
	std::string outputFileName = inputFileName.substr(0,pos) + "_InitialCommunities.txt";
	std::ifstream inputFile(outputFileName.c_str());
	gettimeofday(&startTime,NULL);
	if(inputFile.is_open())
		initialCommunityDefinitionProvidedFileCommunities(g,outputFileName,alpha);
	else
	{
		if(cf == 1.0)
			fuseMethodCalculation(g,outputFileName);
		else
			fuseMethodCalculationMF(g,outputFileName,cf);
	}
	gettimeofday(&endTime,NULL);
	fusingTime = (endTime.tv_sec - startTime.tv_sec)*1000000 + (endTime.tv_usec - startTime.tv_usec);
	//return communitiesBetaFactor;
}
