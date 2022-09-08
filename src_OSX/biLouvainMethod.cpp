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


#include "biLouvainMethod.h"

biLouvainMethod::biLouvainMethod()
{
	_totalPartitioningModularity = 0.0;
	initialCommunityTime = 0.0;
	initialCommunityNeighborsTime = 0.0;
	initialCoClusterMateTime = 0.0;
	candidatesTime = 0.0;
	modularityGainTime = 0.0;
	updateTime = 0.0;
	precalculationCiTime = 0.0;
	precalculationCjTime = 0.0;
	precalculationDTime = 0.0;
	premurataTime = 0.0;
}

biLouvainMethod::~biLouvainMethod(){}

/* AUXILIAR FUNCTIONS AND PROCEDURES */
int biLouvainMethod::findCommunityContainingNode(int nodeId)
{
	int result = -1;
	for(int i=0;i<_numberCommunities;i++)
	{
		for (int j = 0; j < _communities[i].getNumberNodes(); j++)
		{
			if(_communities[i].getNodes()[j]==nodeId)
			{
				result = i;
				break;
			}
		}
	}
	return result;
}

int biLouvainMethod::isRepeated(const std::vector<int> &elements, int newElement)
{
	int result = 0;
	for(unsigned int i =0;i<elements.size();i++)
	{
		if(elements[i]==newElement)
		{
			result = 1;
			break;
		}
	}
	return result;
}

int biLouvainMethod::calculateEdgesBetweenCommunities(Graph &g,int communityCId, int communityDId)
{
	int result = 0;
	for(int i=0;i<_communities[communityCId].getNumberNodes();i++)
	{
		int temp = _communities[communityCId].getNodes()[i];
		for(int j=0;j<_communities[communityDId].getNumberNodes();j++)
		{
			for(int k=0;k<g._graph[_communities[communityDId].getNodes()[j]].getNumberNeighbors();k++)
			{
				if(g._graph[_communities[communityDId].getNodes()[j]].getNeighbors()[k] == temp)
					result++;
			}
		}
	}
	return result;
}

double biLouvainMethod::calculateEdgesBetweenCommunitiesMap(Graph &g,int communityCId, int communityDId)
{
	double result = 0.0;
	for(int i=0;i<_communities[communityCId].getNumberNodes();i++)
		result +=  g._graph[_communities[communityCId].getNodes()[i]].getWeightEdgesToNeighborCommunity(communityDId);
	return result;
}

int biLouvainMethod::findPositionNode(Graph &g,int nodeId)
{
	int result = -1;
	for(int i=0;i<g._numberNodes;i++)
	{
		if(g._graph[i].getId()==nodeId)
		{
			result = i;
			break;
		}
	}
	return result;
}

std::vector<int> biLouvainMethod::findNeighborCommunities(Graph &g,int communityId)
{
	std::vector<int> neighborCommunities;
	int communityContainingNodeId = 0;
	int is_repeated = 0;
	for(int j=0;j<_communities[communityId].getNumberNodes();j++)
	{
		for(int k=0;k<g._graph[_communities[communityId].getNodes()[j]].getNumberNeighbors();k++)
		{
			communityContainingNodeId = findCommunityContainingNode(g._graph[_communities[communityId].getNodes()[j]].getNeighbors()[k]);
			is_repeated = isRepeated(neighborCommunities,communityContainingNodeId);
			if(is_repeated == 0)
				neighborCommunities.push_back(communityContainingNodeId);
		}
	}
	return neighborCommunities;
}

std::vector<int> biLouvainMethod::findNeighborCommunitiesMap(Graph &g,int communityId)
{
	std::vector<int> neighborCommunities;
	std::vector<int> temp;

	for(int i=0;i<_communities[communityId].getNumberNodes();i++)
	{
		temp =  g._graph[_communities[communityId].getNodes()[i]].getNeighborCommunities();
		neighborCommunities.insert(neighborCommunities.end(),temp.begin(),temp.end());
		sort(neighborCommunities.begin(),neighborCommunities.end());
		neighborCommunities.erase(unique( neighborCommunities.begin(), neighborCommunities.end()),neighborCommunities.end());
	}
	return neighborCommunities;
}

std::vector<int> biLouvainMethod::findNeighborCommunitiesWithoutNodeMap(Graph &g,int communityId, int nodeId)
{
	std::vector<int> neighborCommunities;
	std::vector<int> temp;

	for(int i=0;i<_communities[communityId].getNumberNodes()-1;i++)
	{
		temp = g._graph[_communities[communityId].getNodesWithoutNode(nodeId)[i]].getNeighborCommunities();
		neighborCommunities.insert(neighborCommunities.end(), temp.begin(), temp.end() );
		sort(neighborCommunities.begin(),neighborCommunities.end());
		neighborCommunities.erase(unique( neighborCommunities.begin(), neighborCommunities.end()),neighborCommunities.end());
	}
	return neighborCommunities;
}


void biLouvainMethod::calculateCommunitiesModulatiryContribution()
{
	_totalPartitioningModularity = 0.0;
	for(int i=0;i<_numberCommunities;i++)
		_totalPartitioningModularity += _communities[i].getModularityContribution();
	//printf("\n Total partitioning modularity: %f \n", _total_partitioning_modularity);
}

std::vector<int> biLouvainMethod::getDifferentNeighborCommunities(Graph &g,int communityId1, int communityId2)
{
	std::vector<int> result;
	std::vector<int> neighbors_community1;
	std::vector<int> neighbors_community2;
	bool band;
	neighbors_community1 = findNeighborCommunities(g,communityId1);
	neighbors_community2 = findNeighborCommunities(g,communityId2);
	for(unsigned int i=0;i<neighbors_community1.size();i++)
	{
		band = false;
		for(unsigned int j=0;j<neighbors_community2.size();j++)
		{
			if(neighbors_community1[i]==neighbors_community2[j])
			{
				band = true;
				break;
			}
		}
		if(band == false)
			result.push_back(neighbors_community1[i]);
	}
	return result;
}

std::vector<int> biLouvainMethod::getDifferentNeighborCommunitiesMap(Graph &g,int communityId1, int communityId2)
{
	std::vector<int> neighborsCommunity1;
	std::vector<int> neighborsCommunity2;
	std::vector<int>::iterator it;
	neighborsCommunity1 = findNeighborCommunitiesMap(g,communityId1);
	neighborsCommunity2 = findNeighborCommunitiesMap(g,communityId2);
	std::vector<int> result(neighborsCommunity1.size()+neighborsCommunity2.size());
	sort(neighborsCommunity1.begin(),neighborsCommunity1.end());
	sort(neighborsCommunity2.begin(),neighborsCommunity2.end());
	it=set_difference (neighborsCommunity1.begin(),neighborsCommunity1.end(),neighborsCommunity2.begin(),neighborsCommunity2.end(),result.begin());
	result.resize(it-result.begin());
	return result;
}

/* MAIN FUNCTIONS AND PROCEDURES */

//Define initial communities: Each node is a community by itself
void biLouvainMethod::initialCommunityDefinition(Graph &g)
{
	_numberCommunities = g._numberNodes;
	std::unordered_map<int,double> nodesInCommunity;
	for(int i=0;i<g._numberNodes;i++)
	{
		//std::cout<<g._graph[i].getId()<<std::endl;
		g._graph[i].setCommunityId(i);
		nodesInCommunity[g._graph[i].getId()]=g._graph[i].getDegreeNode();
		Community community(i,g._graph[i].getType(),nodesInCommunity,1);
		_communities.push_back(community);
		nodesInCommunity.clear();
	}
}

void biLouvainMethod::initialCommunityDefinitionWithIntraType(Graph &g)
{
        _numberCommunities = g._numberNodes;
        std::unordered_map<int,double> nodesInCommunity;
        std::unordered_map<int,double> nodesIntraType;
	double similarity = 0.0;
        for(int i=0;i<g._numberNodes;i++)
        {
                //std::cout<<g._graph[i].getId()<<std::endl;
                g._graph[i].setCommunityId(i);
                nodesInCommunity[g._graph[i].getId()]=g._graph[i].getDegreeNode();
                nodesIntraType[g._graph[i].getId()]=g._graph[i].getSimilarityNode();//-g._graph[i].getSimilarityIntraTypeNeighbor(i);
                Community community(i,g._graph[i].getType(),nodesInCommunity,nodesIntraType);
                _communities.push_back(community);
		similarity = g._graph[i].getSimilarityIntraTypeNeighbor(i);
                _communities[i].setBetaFactor(calculateCommunityBetaFactor(g,_communities[i].getDescription(),similarity));
		//std::cout<<"\nCommunity: "<<i<<"     Similarity: "<<similarity<<"   Beta Factor AC: "<< _communities[i].getBetaFactor()<<std::endl;
                nodesInCommunity.clear();
                nodesIntraType.clear();
        }
}

//Once a community has been created, it defines the neighbors of that community
void biLouvainMethod::initialCommunityNeighborsDefinition(Graph &g)
{
	int key =0;
	std::unordered_map<int,double> neighborCommunities;
	for(int i=0;i<g._numberNodes;i++)
	{
		for(int j=0;j<g._graph[i].getNumberNeighbors();j++)
		{
			key = g._graph[g._graph[i].getNeighbors()[j]].getCommunityId();
			if(neighborCommunities.find(key)!= neighborCommunities.end())
				neighborCommunities[key] += g._graph[i].getWeightNeighbor(g._graph[g._graph[i].getNeighbors()[j]].getId());
			else
				neighborCommunities[key] = g._graph[i].getWeightNeighbor(g._graph[g._graph[i].getNeighbors()[j]].getId());
		}
		g._graph[i].setNeighborCommunities(neighborCommunities);
		neighborCommunities.clear();
	}
}


//Once a community has been created, it defines the intra type neighbors of that community
void biLouvainMethod::initialIntraTypeCommunityNeighborsDefinition(Graph &g)
{
        int key =0;
        std::unordered_map<int,double> intraTypeNeighborCommunities;
	for(int i=0;i<g._numberNodes;i++)
        {
               	for(int j=0;j<g._graph[i].getNumberIntraTypeNeighbors();j++)
               	{	
                       	key = g._graph[g._graph[i].getIntraTypeNeighbors()[j]].getCommunityId();
	                if(intraTypeNeighborCommunities.find(key)!= intraTypeNeighborCommunities.end())
                                intraTypeNeighborCommunities[key] += g._graph[i].getSimilarityIntraTypeNeighbor(g._graph[g._graph[i].getIntraTypeNeighbors()[j]].getId());
               	        else
                       	        intraTypeNeighborCommunities[key] = g._graph[i].getSimilarityIntraTypeNeighbor(g._graph[g._graph[i].getIntraTypeNeighbors()[j]].getId());
	        }	
                g._graph[i].setIntraTypeNeighborCommunities(intraTypeNeighborCommunities);
               	intraTypeNeighborCommunities.clear();
	}	
}


//If node i has moved to another community we need to update the community to wich i belonged and the community i is moving to.
void biLouvainMethod::updateNodeCommunity(Graph &g,int nodeId, int oldCommunityId, int newCommunityId)
{
	if(_communities[oldCommunityId].getNumberNodes()>0)
	{
		_communities[oldCommunityId].deleteNode(nodeId);
		_communities[newCommunityId].addNode(nodeId,g._graph[nodeId].getDegreeNode());
		g._graph[nodeId].setCommunityId(newCommunityId);
	}
}

void biLouvainMethod::updateNodeIntraTypeCommunity(Graph &g,int nodeId, int oldCommunityId, int newCommunityId)
{
        if(_communities[oldCommunityId].getNumberNodes()>0)
        {
                _communities[oldCommunityId].deleteIntraTypeNode(nodeId);
                _communities[newCommunityId].addIntraTypeNode(nodeId,g._graph[nodeId].getSimilarityNode());
        }
}


//When i moves it has an impact on the communities whose neighbors were CiOld or CiNew. The move could have change the cocluster mate and the contribution
void biLouvainMethod::updateCoClusterMateCommunities(const std::string &changes)
{
	int numberPieces = 0;
	int communityId = 0;
	int pos =  0;
	int aux = 0;
	double contribution = 0.0;
	std::vector<int> coClusterMateCommunities;
	std::string* communities = StringSplitter::split(changes,"$",numberPieces);
	//std::cout << "\n Changes: " << changes << std::endl;
	for(int i = 0;i<numberPieces;i++)
	{
		//printf("%s \n",communities[i].c_str());
		pos =  communities[i].find(":");
		communityId = stoi(communities[i].substr(0,pos));

		//Set the new cocluster mate
		pos =  communities[i].find(",");
		if(pos > -1)
		{
			aux =  communities[i].find(":")+1;
			while(pos != std::string::npos)
			{
				//std::cout << "\nCommunity: " << stoi(communities[i].substr(aux,pos)) << std::endl;
				coClusterMateCommunities.push_back(stoi(communities[i].substr(aux,pos)));
				aux = pos + 1;
				pos =  communities[i].find(",",aux);
			}
			pos =  communities[i].find(":");
			coClusterMateCommunities.push_back(stoi(communities[i].substr(aux,pos)));
		}
		else
		{
			pos =  communities[i].find(":");
			aux =  communities[i].find("#");
			coClusterMateCommunities.push_back(stoi(communities[i].substr(pos+1,aux)));
		}

		_communities[communityId].setCoClusterMateCommunityId(coClusterMateCommunities);
		coClusterMateCommunities.clear();

		//Set the new contribution to modularity
		pos =  communities[i].find("#");
		contribution = stod(communities[i].substr(pos+1,communities[i].length()));
		_communities[communityId].setModularityContribution(contribution);
		//printf("Community ID:%d  Contribution:%f \n" , communityId,contribution);
	}
}

//When i moves it changes the neighbors of CiOld and CiNew as well
void biLouvainMethod::updateNeighborCommunities(Graph &g,int nodeId, int oldCommunityId,int newCommunityId)
{
	for(int i=0;i<g._graph[nodeId].getNumberNeighbors();i++)
	{
		g._graph[g._graph[nodeId].getNeighbors()[i]].deleteNeighborCommunityWeight(oldCommunityId,g._graph[nodeId].getWeightNeighbor(g._graph[nodeId].getNeighbors()[i]));
		g._graph[g._graph[nodeId].getNeighbors()[i]].addNeighborCommunityWeight(newCommunityId,g._graph[nodeId].getWeightNeighbor(g._graph[nodeId].getNeighbors()[i]));
	}
}

//When i moves it changes the beta factor of CiOld and CiNew as well
void biLouvainMethod::updateIntraTypeNeighborCommunities(Graph &g,int nodeId, int oldCommunityId,int newCommunityId)
{
	for(int i=0;i<g._graph[nodeId].getNumberIntraTypeNeighbors();i++)
        {
		g._graph[g._graph[nodeId].getIntraTypeNeighbors()[i]].deleteIntraTypeNeighborCommunitySimilarity(oldCommunityId,g._graph[nodeId].getSimilarityIntraTypeNeighbor(g._graph[nodeId].getIntraTypeNeighbors()[i]));
	        g._graph[g._graph[nodeId].getIntraTypeNeighbors()[i]].addIntraTypeNeighborCommunitySimilarity(newCommunityId,g._graph[nodeId].getSimilarityIntraTypeNeighbor(g._graph[nodeId].getIntraTypeNeighbors()[i]));
	}
}


//Methods for Graph compaction at the end of each Phase. It helps to go out of local maximum

std::unordered_map<int,int> biLouvainMethod::dictionaryCommunitiesNewId()
{
	std::unordered_map<int,int> dictionaryCommunities;
	int id = 0;
	for(int i=0;i<_numberCommunities;i++)
	{
		if(_communities[i].getNumberNodes()>0)
		{
			dictionaryCommunities[_communities[i].getId()] = id;
			 //printf("\nCommunity: %d      Dic: %d", i,id);
			id++;
		}
	}
	return dictionaryCommunities;
}

int biLouvainMethod::numberCommunitiesNonEmpty()
{
	int result = 0;
	for(int i=0;i<_numberCommunities;i++)
	{
		if(_communities[i].getNumberNodes()>0)
			result++;
	}
	return result;
}

std::unordered_map<int,double> biLouvainMethod::compactMetaNodeNeighbors(Graph &g,int &communityId, std::unordered_map<int,int> &dictionaryCommunities)
{
 	std::unordered_map<int,long double> neighborsTemp;
        std::unordered_map<int,long double> errorCalculation;
	std::unordered_map<int,double> neighbors;
 	for(int j=0; j<_communities[communityId].getNumberNodes();j++)
        {
		std::vector<Node> temp = g._graph[_communities[communityId].getNodes()[j]].getNodes();
                nodes.insert(nodes.end(),temp.begin(),temp.end());
                //std::cout << "Community: " << _communities[communityId].getId()<<" Node: "<<_communities[i].getNodes()[j]<<" #Nei: "<<g._graph[_communities[communityId].getNodes()[j]].getNumberNeighbors()<<std::endl;
	        long double y = 0.0, t = 0.0;
        	for(int k=0; k<g._graph[_communities[communityId].getNodes()[j]].getNumberNeighbors();k++)
        	{
        		int idNeighbor = dictionaryCommunities[g._graph[g._graph[_communities[communityId].getNodes()[j]].getNeighbors()[k]].getCommunityId()]; 
                	if(neighborsTemp.find(idNeighbor)!= neighborsTemp.end())
                	{
                		y= g._graph[_communities[communityId].getNodes()[j]].getWeightNeighbor(g._graph[_communities[communityId].getNodes()[j]].getNeighbors()[k])-errorCalculation[idNeighbor];
                        	t = neighborsTemp[idNeighbor] + y;
                        	errorCalculation[idNeighbor] = (t - neighborsTemp[idNeighbor]) - y;
                     	   	neighborsTemp[idNeighbor] = t;
                	}
                	else
                	{
                		neighborsTemp[idNeighbor] = g._graph[_communities[communityId].getNodes()[j]].getWeightNeighbor(g._graph[_communities[communityId].getNodes()[j]].getNeighbors()[k]);
                        	errorCalculation[idNeighbor] = 0.0;
                	}
        	}
	}
       	//errorCalculation.clear();
        //std::stringstream line;
        //line.precision(64);
        //std::cout<<"\nCommunity: " << id;
	for(auto it=neighborsTemp.begin();it!=neighborsTemp.end();++it)
        	neighbors[it->first] = (double)it->second;
        /*      line << "N: " << it->first <<"\t Weight: " << it->second << std::endl;
                std::cout << line.str();
                line.str("");
        }*/
	return neighbors;
}


std::unordered_map<int,double> biLouvainMethod::compactMetaNodeIntraTypeNeighbors(Graph &g,int &communityId, std::unordered_map<int,int> &dictionaryCommunities)
{
        std::unordered_map<int,long double> neighborsTemp;
	std::unordered_map<int,double> neighbors;
        std::unordered_map<int,long double> errorCalculation;
        for(int j=0; j<_communities[communityId].getNumberNodes();j++)
        {
                long double y = 0.0, t = 0.0;
                for(int k=0; k<g._graph[_communities[communityId].getNodes()[j]].getNumberIntraTypeNeighbors();k++)
                {
                        int idNeighbor = dictionaryCommunities[g._graph[g._graph[_communities[communityId].getNodes()[j]].getIntraTypeNeighbors()[k]].getCommunityId()];
                        if(neighborsTemp.find(idNeighbor)!= neighborsTemp.end())
                        {
                                y= g._graph[_communities[communityId].getNodes()[j]].getSimilarityIntraTypeNeighbor(g._graph[_communities[communityId].getNodes()[j]].getIntraTypeNeighbors()[k])-errorCalculation[idNeighbor];
                                t = neighborsTemp[idNeighbor] + y;
                                errorCalculation[idNeighbor] = (t - neighborsTemp[idNeighbor]) - y;
                                neighborsTemp[idNeighbor] = t;
                        }
                        else
                        {
                                neighborsTemp[idNeighbor] = g._graph[_communities[communityId].getNodes()[j]].getSimilarityIntraTypeNeighbor(g._graph[_communities[communityId].getNodes()[j]].getIntraTypeNeighbors()[k]);
                                errorCalculation[idNeighbor] = 0.0;
                        }
                }
        }
	for(auto it=neighborsTemp.begin();it!=neighborsTemp.end();++it)
        	neighbors[it->first] = (double)it->second;
        return neighbors;
}

void biLouvainMethod::fromCommunitiesToNodes(Graph &g)
{
	int lastIdPartitionV1 = -1;
	std::unordered_map<int,int> dictionaryCommunities = dictionaryCommunitiesNewId();
	int numberNodes = dictionaryCommunities.size();
	MetaNode*_newGraph = new MetaNode[numberNodes];
	std::unordered_map<int,double> neighbors;
	//std::cout << "\nPrint del from: "<< numberNodes << "\t " << dictionaryCommunities.size() << std::endl;
	int id =0;
	for(int i=0;i<_numberCommunities;i++)
	{
		if(_communities[i].getNumberNodes()>0)
		{
			neighbors = compactMetaNodeNeighbors(g,i,dictionaryCommunities);
			MetaNode metanode(id,_communities[i].getDescription(),nodes,neighbors,-1);
			_newGraph[id] = metanode;
			if(_alpha != 1.0)
			{
				neighbors.clear();
				neighbors = compactMetaNodeIntraTypeNeighbors(g,i,dictionaryCommunities);
				_newGraph[id].setIntraTypeNeighbors(neighbors);
			//	printf("\nPrint del from:::Community: %d    BF: %f ",i,_communities[i].getBetaFactor());
			}
			if(metanode.getType()=="V1")
				lastIdPartitionV1++;
			id++;
			nodes.clear();
			neighbors.clear();
		}
	}
	Graph compactedGraph(_newGraph,numberNodes,g._numberEdges,g._weightEdges,g._weightEdgesV1,g._weightEdgesV2,lastIdPartitionV1);
	if(_alpha != 1.0)
	{
		compactedGraph.setLambdaV1(g._lambdaV1);
		compactedGraph.setLambdaV2(g._lambdaV2);
		compactedGraph.setSimilarityV1(g._sumSimilarityV1);
		compactedGraph.setSimilarityV2(g._sumSimilarityV2);
	}
	g.destroyGraph();
	g = compactedGraph;
	dictionaryCommunities.clear();
	_communities.clear();
	//for(int i=0;i<g._numberNodes;i++)
	//	std::cout << g._graph[i].getId() << "  " << g._graph[i].getNumberNodes() << std::endl;
}


//MAIN METHOD FOR BILOUVAIN ALGORITHM

//Returns an array hat contains the order in which the nodes should be analyzed during each iteration 
int* biLouvainMethod::nodesOrderToProcess(Graph &g,int optionOrder)
{
	int* nodesOrder = new int[g._numberNodes];
	switch(optionOrder){
		case 1:
		{
			for(int i=0;i<g._numberNodes;i++)
				nodesOrder[i] = i;
			break;
		}
		case 2:
		{
			int k = 1;
			int l = 1;
			int m = 1;
			int numberNodesV2 = g._numberNodes - (g._lastIdPartitionV1+1);
			for(int i=0;i<g._numberNodes;i++)   //Creation of the nodes
			{
				if((g._lastIdPartitionV1 + 1) ==  numberNodesV2)
				{	if(i<=g._lastIdPartitionV1)
						nodesOrder[i*2] = i;
					else
					{
						nodesOrder[(i-(g._lastIdPartitionV1+1))+k]=i;
						k++;
					}
				}
				else if((g._lastIdPartitionV1 + 1) >  numberNodesV2)
				{
					if(i<=g._lastIdPartitionV1)
					{
						if(k<=numberNodesV2+1)
						{	nodesOrder[i*2] = i;
							m = (k*2)-1;
							k++;
						}
						else
						{
							nodesOrder[m]=i;
							m++;
						}
					}else
					{
						nodesOrder[(i-(g._lastIdPartitionV1+1))+l]=i;
						l++;
					}
				}
				else if((g._lastIdPartitionV1 + 1) <  numberNodesV2)
				{
					if(i<=g._lastIdPartitionV1)
						nodesOrder[i*2]=i;
					else
					{
						if(l<=g._lastIdPartitionV1 + 1)
						{
							nodesOrder[(i-(g._lastIdPartitionV1+1))+l]=i;
							m=(l*2);
							l++;
						}
						else
						{
							nodesOrder[m]=i;
							m++;
						}
					}
				}
			}
			break;
		}
		case 3:
		{
			//srand(time(NULL));
			int positions[g._numberNodes];
			for(int i=0;i<g._numberNodes;i++)   //Initialize array elements to 0
				positions[i] = 0;
			for(int i=0;i<g._numberNodes;i++)   //Select the order randomly with the same seed for rand
			{
				bool band = true;
				int number = 0;
				while(band)
				{
					number = rand() % g._numberNodes;
					if(positions[number] == 0)
					{
						positions[number] = 1;
						band = false;
						nodesOrder[i] = number;
					}
				}
			}
			break;
		}
	}
	return nodesOrder;
}


newDataCommunity biLouvainMethod::calculateDeltaGainModularity(Graph &g,MetaNode &node,int &communityId, int newCommunityId,int option)
{
        double deltaModularityContribution = 0.0;
        newDataCommunity newCalculationModulatiry;

        newCalculationModulatiry = CoClusterMateDefinitionPrecalculation(g,node,communityId, newCommunityId,option);
        deltaModularityContribution = newCalculationModulatiry.newModularityContribution - _communities[communityId].getModularityContribution();
        //printf("Before: %f  After:%f  Difference:%f \n",_communities[communityId].getModularityContribution(),newCalculationModulatiry.newModularityContribution,deltaModularityContribution);
        newCalculationModulatiry.newModularityContribution = deltaModularityContribution;
        return newCalculationModulatiry;
}

double biLouvainMethod::calculateMaxModularityGainIteration(Graph &g,int* &nodesOrderExecution)
{
	struct timeval t1,t2,t3,t4,t5,t6;
	double maxModularityGainIteration = 0.0;  
	std::vector<int> candidates;
	double lambda = 0.0;
	double lambdaD = 0.0;
	double betaFactorCurrentCommunity = 0.0;
	double betaFactorCandidateCommunity = 0.0;
	for(int i=0;i<g._numberNodes;i++)
	{
		std::stringstream gainString;
		double gainDoble = 0.0;
		std::stringstream changes;
		newDataCommunity deltaModularityGain;
		double maxDeltaModularityGain = -1.0;
		double totalDeltaModularityGain = 0.0;
		int currentCommunity = g._graph[nodesOrderExecution[i]].getCommunityId();
		if(_alpha != 1.0)
		{
			if(_communities[currentCommunity].getDescription()=="V1")
			{
				lambda = g._lambdaV1;
				lambdaD = g._lambdaV2;
			}
			else
			{
				lambda = g._lambdaV2;
				lambdaD = g._lambdaV1;
			}
		}
		//Find CANDIDATE COMMUNITIES to which node i can move to
		gettimeofday(&t1,NULL);
		std::vector<int> temp;
		if(_alpha != 0.0)
		{
			for(int j=0;j<g._graph[nodesOrderExecution[i]].getNumberNeighbors();j++)
			{
				temp = g._graph[g._graph[nodesOrderExecution[i]].getNeighbors()[j]].getNeighborCommunities();
				candidates.insert(candidates.end(),temp.begin(),temp.end());
				temp.clear();
			}
		}
		if(_alpha != 1.0) 
		{
			temp = g._graph[nodesOrderExecution[i]].getIntraTypeNeighborCommunities();
                        candidates.insert(candidates.end(),temp.begin(),temp.end());
                        temp.clear();
		}	
		sort(candidates.begin(),candidates.end());
		candidates.erase(unique(candidates.begin(),candidates.end()),candidates.end());
		std::vector<int>::iterator position = find(candidates.begin(),candidates.end(),currentCommunity);
		if (position != candidates.end()) candidates.erase(position);
		gettimeofday(&t2,NULL);

		/*printf("\n Node: %d \n",_graph[nodesOrderExecution[i]].getId());
		StringSplitter::printVector(candidates);
		printf("\n");*/

		//Calculate the GAIN IN MODULARITY of the new setup of the structure if node i actually moves
		gettimeofday(&t3,NULL);
		changes.str("");
		//changes.precision(64);
		//printf("1 \n");
                //Calculate Delta QB for the community to which i belongs to (Ci)
		deltaModularityGain = calculateDeltaGainModularity(g,g._graph[nodesOrderExecution[i]],currentCommunity,0,1);
		betaFactorCurrentCommunity = deltaModularityGain.newBetaFactor;
		gainDoble = deltaModularityGain.newModularityContribution + _communities[currentCommunity].getModularityContribution();
		changes << deltaModularityGain.coClusterMateCommunityId << "#" << gainDoble << "$";
		//gainString << gainDoble;
		//totalDeltaModularityGain += deltaModularityGain.newModularityContribution - (gainDoble - stof(gainString.str()));
		totalDeltaModularityGain += deltaModularityGain.newModularityContribution;
		//Calculate Delta QB for the set of candidate communities (Cj)
		int candidateCommunity = -1;
		std::stringstream maxChangesCandidate;
		if(candidates.size()>0)
		{
			std::stringstream temp;
			double betaF = 0.0;
			double candidateDeltaModularityGain = 0.0;
			//betaFactorCandidateCommunity = 0.0;
			betaFactorCandidateCommunity = _communities[currentCommunity].getBetaFactor();
			for(unsigned int j=0;j<candidates.size();j++)
			{
				candidateDeltaModularityGain = 0.0;				
				betaF = 0.0;
				temp.str("");
				//temp.precision(64);
				//printf("2 \n"); 
				deltaModularityGain = calculateDeltaGainModularity(g,g._graph[nodesOrderExecution[i]],candidates[j],currentCommunity,2);
				betaF = deltaModularityGain.newBetaFactor;
				gainDoble = deltaModularityGain.newModularityContribution + _communities[candidates[j]].getModularityContribution();
				temp << deltaModularityGain.coClusterMateCommunityId << "#" << gainDoble << "$";
                                //gainString << gainDoble;
				//candidateDeltaModularityGain += deltaModularityGain.newModularityContribution-(gainDoble - stof(gainString.str()));
				 candidateDeltaModularityGain += deltaModularityGain.newModularityContribution;
				//Calculate Delta QB for neighbors of candidate communities (Dj)
				//printf("4 \n");
				if((betaF >= lambda)&&(betaF > betaFactorCandidateCommunity))
				{
					std::vector<int>neighborCommunities = findNeighborCommunitiesMap(g,candidates[j]);
					for(unsigned int k=0;k<neighborCommunities.size();k++)
					{	
						deltaModularityGain = calculateDeltaGainModularity(g,g._graph[nodesOrderExecution[i]],neighborCommunities[k],candidates[j],4);
						gainDoble = deltaModularityGain.newModularityContribution + _communities[neighborCommunities[k]].getModularityContribution();
						temp << deltaModularityGain.coClusterMateCommunityId << "#" << gainDoble << "$";
						//gainString.str("");
						//gainString << gainDoble;
					        //candidateDeltaModularityGain += deltaModularityGain.newModularityContribution-(gainDoble - stof(gainString.str()));
						candidateDeltaModularityGain += deltaModularityGain.newModularityContribution;
					}
					if((candidateDeltaModularityGain > maxDeltaModularityGain)&&(betaF > _communities[currentCommunity].getBetaFactor()))
					{
						maxDeltaModularityGain = candidateDeltaModularityGain;
						candidateCommunity = candidates[j];
						maxChangesCandidate.str("");
						maxChangesCandidate << temp.str();
						if(_alpha != 1.0)
							betaFactorCandidateCommunity = betaF;
					}
				}
			}
			candidates.clear();
			changes << maxChangesCandidate.str();
			maxChangesCandidate.str("");
			totalDeltaModularityGain += maxDeltaModularityGain;
		}
		if((candidateCommunity != -1)&&(betaFactorCandidateCommunity >= lambda))
		{	
			//Calculate Delta QB for the neighbors of Ci (Di)		
			std::vector<int> differentNeighborCommunities = getDifferentNeighborCommunitiesMap(g,currentCommunity,candidateCommunity);
			if(differentNeighborCommunities.size()>0)
			{
				//printf("3 \n");
				for(unsigned int j=0;j<differentNeighborCommunities.size();j++)
				{
					deltaModularityGain = calculateDeltaGainModularity(g,g._graph[nodesOrderExecution[i]],differentNeighborCommunities[j],candidateCommunity,3);
					//printf("  Delta Modularity Gain: %f \n", deltaModularityGain);
					gainDoble = deltaModularityGain.newModularityContribution + _communities[differentNeighborCommunities[j]].getModularityContribution();
					changes << deltaModularityGain.coClusterMateCommunityId << "#" << gainDoble << "$";
	                                //gainString << gainDoble;
	        	       	        //totalDeltaModularityGain += deltaModularityGain.newModularityContribution - (gainDoble- stof(gainString.str()));
					totalDeltaModularityGain += deltaModularityGain.newModularityContribution;
				}
			}
			changes.str(changes.str().substr(0,changes.str().length()-1));
			std::string changesMade = changes.str();
			//printf("Delta Modularity Gain: %f  Changes: %s \n", totalDeltaModularityGain, changesMade.c_str());
			gettimeofday(&t4,NULL);

			//If the Gain in modularity
			gettimeofday(&t5,NULL);
			if(totalDeltaModularityGain > 0.0)
			{
				updateNodeCommunity(g,g._graph[nodesOrderExecution[i]].getId(),currentCommunity,candidateCommunity);
				updateCoClusterMateCommunities(changesMade);
				updateNeighborCommunities(g,g._graph[nodesOrderExecution[i]].getId(),currentCommunity,candidateCommunity);
				if(_alpha != 1.0)
				{
				    _communities[candidateCommunity].setBetaFactor(betaFactorCandidateCommunity);
				    _communities[currentCommunity].setBetaFactor(betaFactorCurrentCommunity);
				    updateNodeIntraTypeCommunity(g,g._graph[nodesOrderExecution[i]].getId(),currentCommunity,candidateCommunity);
				    updateIntraTypeNeighborCommunities(g,g._graph[nodesOrderExecution[i]].getId(),currentCommunity,candidateCommunity);
				}
				printf("\n Node: %d  From Community: %d  To Community: %d  Maximum Modularity Gain: %.15lf",g._graph[nodesOrderExecution[i]].getId(), currentCommunity,candidateCommunity,totalDeltaModularityGain);
				maxModularityGainIteration += totalDeltaModularityGain;
				//printCommunitiesContributionModularity();
			}
			gettimeofday(&t6,NULL);
			candidatesTime += (t2.tv_sec - t1.tv_sec)*1000000 + (t2.tv_usec - t1.tv_usec);
			modularityGainTime += (t4.tv_sec - t3.tv_sec)*1000000 + (t4.tv_usec - t3.tv_usec);
			updateTime += (t6.tv_sec - t5.tv_sec)*1000000 + (t6.tv_usec - t5.tv_usec);
		}
		double nodeTime = (t6.tv_sec - t1.tv_sec)*1000000 + (t6.tv_usec - t1.tv_usec);
	   //}
	}
	return maxModularityGainIteration;
}




double biLouvainMethod::calculateMaxModularityGainIterationIntraType(Graph &g,int* &nodesOrderExecution)
{
        struct timeval t1,t2,t3,t4,t5,t6;
        double maxModularityGainIteration = 0.0;
        std::vector<int> candidates;
        double lambda = 0.0;
        double lambdaD = 0.0;
        double betaFactorCurrentCommunity = 0.0;
        double betaFactorCandidateCommunity = 0.0;
	double newContributionCurrentCommunity = 0.0;
	double newContributionCandidateCommunity = 0.0;
        for(int i=0;i<g._numberNodes;i++)
        {
                double gainDoble = 0.0;
                newDataCommunity deltaModularityGain;
                double maxDeltaModularityGain = -1.0;
                double totalDeltaModularityGain = 0.0;
                int currentCommunity = g._graph[nodesOrderExecution[i]].getCommunityId();
                if(_communities[currentCommunity].getDescription()=="V1")
                {
                	lambda = g._lambdaV1;
                        lambdaD = g._lambdaV2;
                }
                else
                {
                        lambda = g._lambdaV2;
                        lambdaD = g._lambdaV1;
                }
                //Find CANDIDATE COMMUNITIES to which node i can move to
                gettimeofday(&t1,NULL);
                candidates = g._graph[nodesOrderExecution[i]].getIntraTypeNeighborCommunities();
                sort(candidates.begin(),candidates.end());
                std::vector<int>::iterator position = find(candidates.begin(),candidates.end(),currentCommunity);
                if (position != candidates.end()) candidates.erase(position);
                gettimeofday(&t2,NULL);
                //Calculate the GAIN IN MODULARITY of the new setup of the structure if node i actually moves
                gettimeofday(&t3,NULL);
                //Calculate Delta QB for the community to which i belongs to (Ci)
                deltaModularityGain = calculateDeltaGainModularity(g,g._graph[nodesOrderExecution[i]],currentCommunity,0,1);
                betaFactorCurrentCommunity = deltaModularityGain.newBetaFactor;
                newContributionCurrentCommunity = deltaModularityGain.newModularityContribution + _communities[currentCommunity].getModularityContribution();
                totalDeltaModularityGain += deltaModularityGain.newModularityContribution;
		//Calculate Delta QB for the set of candidate communities (Cj)
                int candidateCommunity = -1;
                if(candidates.size()>0)
                {
                        double betaF = 0.0;
			double candidateDeltaModularityGain = 0.0;
			betaFactorCandidateCommunity = _communities[currentCommunity].getBetaFactor();
                        for(unsigned int j=0;j<candidates.size();j++)
                        {
                                candidateDeltaModularityGain = 0.0;
				betaF = 0.0;
                                deltaModularityGain = calculateDeltaGainModularity(g,g._graph[nodesOrderExecution[i]],candidates[j],currentCommunity,2);
                                betaF = deltaModularityGain.newBetaFactor;
                                gainDoble = deltaModularityGain.newModularityContribution + _communities[candidates[j]].getModularityContribution();
                                candidateDeltaModularityGain += deltaModularityGain.newModularityContribution;
                                if((betaF >= lambda)&&(betaF > betaFactorCandidateCommunity))
                                {
                                	maxDeltaModularityGain = candidateDeltaModularityGain;
                                        candidateCommunity = candidates[j];
                                        betaFactorCandidateCommunity = betaF;
					newContributionCandidateCommunity = gainDoble;
                                }
                        }
			candidates.clear();
                        totalDeltaModularityGain += maxDeltaModularityGain;
                }
                if((candidateCommunity != -1)&&(betaFactorCandidateCommunity >= lambda))
                {
			gettimeofday(&t4,NULL);
                        //If the Gain in modularity
                        gettimeofday(&t5,NULL);
                        if(totalDeltaModularityGain > 0.0)
                        {
                             updateNodeCommunity(g,g._graph[nodesOrderExecution[i]].getId(),currentCommunity,candidateCommunity);
                             updateNeighborCommunities(g,g._graph[nodesOrderExecution[i]].getId(),currentCommunity,candidateCommunity);
			     _communities[candidateCommunity].setBetaFactor(betaFactorCandidateCommunity);
                             _communities[currentCommunity].setBetaFactor(betaFactorCurrentCommunity);
			     _communities[candidateCommunity].setModularityContribution(newContributionCandidateCommunity);
                             _communities[currentCommunity].setModularityContribution(newContributionCurrentCommunity);	
                             updateNodeIntraTypeCommunity(g,g._graph[nodesOrderExecution[i]].getId(),currentCommunity,candidateCommunity);
                             updateIntraTypeNeighborCommunities(g,g._graph[nodesOrderExecution[i]].getId(),currentCommunity,candidateCommunity);
                             printf("\n Node: %d  From Community: %d  To Community: %d  Maximum Modularity Gain: %.15lf",g._graph[nodesOrderExecution[i]].getId(), currentCommunity,candidateCommunity,totalDeltaModularityGain);
			     maxModularityGainIteration += totalDeltaModularityGain;
                             //printCommunitiesContributionModularity();
                        }
                        gettimeofday(&t6,NULL);
                        candidatesTime += (t2.tv_sec - t1.tv_sec)*1000000 + (t2.tv_usec - t1.tv_usec);
                        modularityGainTime += (t4.tv_sec - t3.tv_sec)*1000000 + (t4.tv_usec - t3.tv_usec);
                        updateTime += (t6.tv_sec - t5.tv_sec)*1000000 + (t6.tv_usec - t5.tv_usec);
                }
                double nodeTime = (t6.tv_sec - t1.tv_sec)*1000000 + (t6.tv_usec - t1.tv_usec);
           //}
        }
        return maxModularityGainIteration;
}




void biLouvainMethod::biLouvainMethodAlgorithm(Graph &g,double cutoffIterations, double cutoffPhase, int optionOrder,std::unordered_map<int,std::string> &bipartiteOriginalEntities,const std::string &inputFileName,const std::string &outputFileName,double &alpha)
{
	_alpha = alpha;
	//Create output files: one for the modularity gain during each iteration and one for the obtained communities after the iterations are completed
	int pos=0;
	if(outputFileName.empty())
	{
	     pos = inputFileName.find_last_of(".");
	     _outputFileName = inputFileName.substr(0,pos);
	}
	else
	{
		pos = outputFileName.find_last_of(".");
		_outputFileName = outputFileName.substr(0,pos);
	}	
        std::string outputModularityGain = _outputFileName + "_ResultsModularity.txt";
	std::ofstream outfileMG;
	outfileMG.open(outputModularityGain.c_str(),std::ios::out|std::ios::trunc);

	initialCommunityTime = 0.0;
	initialCommunityNeighborsTime = 0.0;
	initialCoClusterMateTime = 0.0;
	candidatesTime = 0.0;
	modularityGainTime = 0.0;
	updateTime = 0.0;
	precalculationCiTime = 0.0;
	precalculationCjTime = 0.0;
	precalculationDTime = 0.0;
	struct timeval t7,t8,t9,t10,t11,t12;
	std::stringstream line;
	double totalModularity = 0.0;
	double phaseModularity = 1;
	int phases = 1;
	int* nodesOrderExecution = NULL;

	//PHASE
	while((phaseModularity-totalModularity) > cutoffPhase)
	{
		//INITIALIZATION STEPS
		gettimeofday(&t7,NULL);
		if(alpha != 1.0)
			initialCommunityDefinitionWithIntraType(g);
		else
			initialCommunityDefinition(g);
		gettimeofday(&t8,NULL);
		initialCommunityTime += (t8.tv_sec - t7.tv_sec)*1000000 + (t8.tv_usec - t7.tv_usec);
		gettimeofday(&t9,NULL);
		initialCommunityNeighborsDefinition(g);
		if(alpha != 1.0)
			initialIntraTypeCommunityNeighborsDefinition(g);
		gettimeofday(&t10,NULL);
		initialCommunityNeighborsTime += (t10.tv_sec - t9.tv_sec)*1000000 + (t10.tv_usec - t9.tv_usec);
		gettimeofday(&t11,NULL);
		totalModularity = CoClusterMateDefinitionAllCommunities(g,0,_numberCommunities);
		phaseModularity = totalModularity;
		gettimeofday(&t12,NULL);
		initialCoClusterMateTime += (t12.tv_sec - t11.tv_sec)*1000000 + (t12.tv_usec - t11.tv_usec);

		printf("\n\n ::: Phase %d :::", phases);
		line.str("");
		line << "--- Phase: " <<  phases << "\n";
		outfileMG << line.str();
		printf("\n Initial Total partitioning modularity: %.15lf", totalModularity);
		line.str("");
		line.precision(15);
		line << "Initial Total Modularity: " <<  totalModularity << "\n";
		outfileMG << line.str();
		printf("\n\n ::: Initial Communities :::");
		printCommunities(g);
		//printCommunitiesContributionModularity();
		int iterations = 1;
		double _cutoffIterations = 2.0;

		nodesOrderExecution = nodesOrderToProcess(g,optionOrder);
		//for(int i=0;i<g._numberNodes;i++)
		//{
		//	printf("%d \t %d \n",i,nodesOrderExecution[i]);
		//}

		//ITERATION		
		while(_cutoffIterations > cutoffIterations)
		{
			printf("\n\n ::: Iteration: %d Start :::",iterations);
			double maxModularityGainIteration = calculateMaxModularityGainIteration(g,nodesOrderExecution);
			calculateCommunitiesModulatiryContribution();
			printf("\n\n ::: Iteration: %d End  :::  Maximum Modularity Gain: %.15lf", iterations,maxModularityGainIteration);
			line.str("");
			line.precision(15);
			line << "Iteration: " <<  iterations << " - Maximum Modularity Gain: " << maxModularityGainIteration << "\n";
			outfileMG << line.str();
			_cutoffIterations = maxModularityGainIteration;
			phaseModularity += maxModularityGainIteration;
			iterations++;
		}
		if((phaseModularity-totalModularity) > cutoffPhase)
		{
			//printCommunitiesContributionModularity();
			fromCommunitiesToNodes(g);
			//printCommunitiesContributionModularity();
			phases++;
		}
	}	
	line.str("");
	line.precision(15);
	line << "\n--- Final Murata+ Modularity: " <<  totalModularity;
	outfileMG <<line.str();
	outfileMG.close();
	printAllCommunityNodeswithSingletons(g,bipartiteOriginalEntities);
	printCoClusterCommunitiesFile();
	_communities.clear();
	delete[] nodesOrderExecution;
}



void biLouvainMethod::biLouvainMethodAlgorithmIntraType(Graph &g,double cutoffIterations, double cutoffPhase, int optionOrder,std::unordered_map<int,std::string> &bipartiteOriginalEntities,const std::string &inputFileName,const std::string &outputFileName)
{
        _alpha = 0.0;
        //Create output files: one for the modularity gain during each iteration and one for the obtained communities after the iterations are completed
        int pos=0;
        if(outputFileName.empty())
        {
             pos = inputFileName.find_last_of(".");
             _outputFileName = inputFileName.substr(0,pos);
        }
        else
        {
                pos = outputFileName.find_last_of(".");
                _outputFileName = outputFileName.substr(0,pos);
        }
        std::string outputModularityGain = _outputFileName + "_ResultsModularity.txt";
        std::ofstream outfileMG;
        outfileMG.open(outputModularityGain.c_str(),std::ios::out|std::ios::trunc);

        initialCommunityTime = 0.0;
        initialCommunityNeighborsTime = 0.0;
        initialCoClusterMateTime = 0.0;
        candidatesTime = 0.0;
        modularityGainTime = 0.0;
        updateTime = 0.0;
        precalculationCiTime = 0.0;
        precalculationCjTime = 0.0;
        precalculationDTime = 0.0;
        struct timeval t7,t8,t9,t10,t11,t12;
        std::stringstream line;
        double totalModularity = 0.0;
        double phaseModularity = 1;
        int phases = 1;
        int* nodesOrderExecution = NULL;

        //PHASE
        while((phaseModularity-totalModularity) > cutoffPhase)
        {
                //INITIALIZATION STEPS
                gettimeofday(&t7,NULL);
                initialCommunityDefinitionWithIntraType(g);
                gettimeofday(&t8,NULL);
                initialCommunityTime += (t8.tv_sec - t7.tv_sec)*1000000 + (t8.tv_usec - t7.tv_usec);
                gettimeofday(&t9,NULL);
                initialCommunityNeighborsDefinition(g);
                initialIntraTypeCommunityNeighborsDefinition(g);
                gettimeofday(&t10,NULL);
                initialCommunityNeighborsTime += (t10.tv_sec - t9.tv_sec)*1000000 + (t10.tv_usec - t9.tv_usec);
                gettimeofday(&t11,NULL);
                totalModularity = IntraTypeDefinitionAllCommunities(g,0,_numberCommunities);
                phaseModularity = totalModularity;
                gettimeofday(&t12,NULL);
                initialCoClusterMateTime += (t12.tv_sec - t11.tv_sec)*1000000 + (t12.tv_usec - t11.tv_usec);

                printf("\n\n ::: Phase %d :::", phases);
                line.str("");
                line << "--- Phase: " <<  phases << "\n";
                outfileMG << line.str();
                printf("\n Initial Total partitioning modularity: %.15lf", totalModularity);
                line.str("");
                line.precision(15);
                line << "Initial Total Modularity: " <<  totalModularity << "\n";
                outfileMG << line.str();
                printf("\n\n ::: Initial Communities :::");
                //printCommunities(g);
                //printCommunitiesContributionModularity();
                int iterations = 1;
                double _cutoffIterations = 2.0;

                nodesOrderExecution = nodesOrderToProcess(g,optionOrder);
                //for(int i=0;i<g._numberNodes;i++)
                //{
                //      printf("%d \t %d \n",i,nodesOrderExecution[i]);
                //}

                //ITERATION             
                while(_cutoffIterations > cutoffIterations)
                {
                        printf("\n\n ::: Iteration: %d Start :::",iterations);
                        double maxModularityGainIteration = calculateMaxModularityGainIterationIntraType(g,nodesOrderExecution);
                        calculateCommunitiesModulatiryContribution();
                        printf("\n\n ::: Iteration: %d End  :::  Maximum Modularity Gain: %.15lf", iterations,maxModularityGainIteration);
                        line.str("");
                        line.precision(15);
                        line << "Iteration: " <<  iterations << " - Maximum Modularity Gain: " << maxModularityGainIteration << "\n";
                        outfileMG << line.str();
                        _cutoffIterations = maxModularityGainIteration;
                        phaseModularity += maxModularityGainIteration;
                        iterations++;
                }
                if((phaseModularity-totalModularity) > cutoffPhase)
                {
                        //printCommunitiesContributionModularity();
                        fromCommunitiesToNodes(g);
		        phases++;
                }
        }
        line.str("");
        line.precision(15);
        line << "\n--- Final Murata+ Modularity: " <<  totalModularity;
        outfileMG <<line.str();
        outfileMG.close();
        printAllCommunityNodeswithSingletons(g,bipartiteOriginalEntities);
        _communities.clear();
        delete[] nodesOrderExecution;
}

//Showing information

int biLouvainMethod::numberNodesInsideCommunity(Graph &g,int communityId)
{
	int result = 0;
	if(_communities[communityId].getNumberNodes()>0)
	{
		for(int i=0; i<_communities[communityId].getNumberNodes();i++)
			result += g._graph[_communities[communityId].getNodes()[i]].getNumberNodes();
	}
	return result;
}

void biLouvainMethod::generateOutputFile(std::string text, std::string fileName)
{
	std::ofstream file;
	file.open(fileName.c_str(), std::ios::out|std::ios::app);
	file << text;
	file.close();
}


void biLouvainMethod::printCommunitiesContributionModularity()
{
	for(int i=0;i<_numberCommunities;i++)
		printf("\nCommunity ID: %d  Contribution: %f  BF: %f\n", _communities[i].getId(),_communities[i].getModularityContribution(),_communities[i].getBetaFactor());
}

void biLouvainMethod::printCommunities(Graph &g)
{
	printf("\n");
	for(int i=0;i<_numberCommunities;i++)
	{
		printf("\nCommunity ID: %d #Nodes: %d  Correspondent Community ID:", _communities[i].getId(),numberNodesInsideCommunity(g,_communities[i].getId()));
		for(unsigned int j=0;j<_communities[i].getCoClusterMateCommunityId().size();j++)
			printf("  %d",_communities[i].getCoClusterMateCommunityId()[j]);
	}
}


void biLouvainMethod::printCoClusterCommunitiesFile()
{
	int countCoClusters = 0;
	std::stringstream line;
	std::string outputCoclusters = _outputFileName + "_ResultsCoClusterCommunities.txt";
	std::ofstream outfileCCC;
	outfileCCC.open(outputCoclusters.c_str(),std::ios::out|std::ios::trunc);
	for(int i=0;i<_numberCommunities;i++)
	{
		countCoClusters++;
		line.str("");
		//line <<"\nCommunity ID: " << _communities[i].getId() << "\t#Nodes: " << numberNodesInsideCommunity(_communities[i].getId()) << "\tCorrespondent Community ID: ";
		line << "\nCoCluster " << countCoClusters << ":"<< _communities[i].getDescription() << "(" << _communities[i].getId() << ")-";
		for(unsigned int j=0;j<_communities[i].getCoClusterMateCommunityId().size();j++)
			line << _communities[i].getCoClusterMateCommunityId()[j] << "  ";
		outfileCCC << line.str().substr(0,line.str().length()-2);
	}
	outfileCCC.close();
}

void biLouvainMethod::printAllCommunityNodes(Graph &g)
{
	std::string outputCommunities = _outputFileName + "_ResultsCommunities.txt";
	std::ofstream outfileC;
	outfileC.open(outputCommunities.c_str(),std::ios::out|std::ios::trunc);
	int singletonsV1 = 0 ;
	int singletonsV2 = 0 ;
	int cont = 0;
	std::stringstream line;
	for(int i=0;i<_numberCommunities;i++)
	{
		if (!_communities[i].getNodes().empty())
		{
			if(_communities[i].getNumberNodes()==1)
			{
				if(g._graph[_communities[i].getNodes()[0]].getNumberNodes()==1)
				{
					if(_communities[i].getDescription() == "V1") singletonsV1++;
					else singletonsV2++;
					std::cout << "Singleton: " << g._graph[_communities[i].getNodes()[0]].getNodes()[0].getIdInput() << std::endl;
				}
				else
				{
					line.str("");
					line << "Community " << cont++ << "[" << _communities[i].getDescription() << "]: ";
					outfileC << line.str();
					line.str("");
					for(int k=0;k<g._graph[_communities[i].getNodes()[0]].getNumberNodes();k++)
						line << g._graph[_communities[i].getNodes()[0]].getNodesSorted()[k].getIdInput() << ",";
					outfileC << line.str().substr(0,line.str().length()-1) << "\n";
				}
			}
			else
			{
				line.str("");
				line << "Community " << cont++ << "[" << _communities[i].getDescription() << "]: ";
				outfileC << line.str();
				line.str("");
				for(int j=0;j<_communities[i].getNumberNodes();j++)
				{
					for(int k=0;k< g._graph[_communities[i].getNodes()[j]].getNumberNodes();k++)
						line << g._graph[_communities[i].getNodes()[j]].getNodesSorted()[k].getIdInput() << ",";
				}
				outfileC << line.str().substr(0,line.str().length()-1) << "\n";
			}
		}
	}
	line.str("");
	line << "\nSingletons Partition V1: " << singletonsV1;
	outfileC << line.str();
	line.str("");
	line << "\nSingletons Partition V2: " << singletonsV2;
	outfileC << line.str();
	outfileC.close();
}

std::string biLouvainMethod::listNodesCommunities(Graph &g)
{
	std::stringstream line;
	std::unordered_map<int,int> nodes;
	line.str("");
	for(int i=0;i<_numberCommunities;i++)
	{
		for(int j=0;j<_communities[i].getNumberNodes();j++)
		{
			for(int k=0;k< g._graph[_communities[i].getNodes()[j]].getNumberNodes();k++)
				nodes[g._graph[_communities[i].getNodes()[j]].getNodesSorted()[k].getIdInput()]=i;
		}
	}
	std::map<int,int> ordered(nodes.begin(),nodes.end());
	for(auto it=ordered.begin();it!=ordered.end();++it)
		line << it->second << ",";
	return line.str().substr(0,line.str().length()-1);
}

void biLouvainMethod::printAllCommunityNodeswithSingletons(Graph &g,std::unordered_map<int,std::string> &bipartiteOriginalEntities)
{
	std::string outputCommunities = _outputFileName + "_ResultsCommunities.txt";
	std::ofstream outfileC;
	outfileC.open(outputCommunities.c_str(),std::ios::out|std::ios::trunc);
	int singletonsV1 = 0 ;
	int singletonsV2 = 0 ;
	int cont = 0;
	std::stringstream line;
	for(int i=0;i<_numberCommunities;i++)
	{
		if (!_communities[i].getNodes().empty())
		{
			if(_communities[i].getNumberNodes()==1)
			{
				if(g._graph[_communities[i].getNodes()[0]].getNumberNodes()==1)
				{
					if(_communities[i].getDescription() == "V1") singletonsV1++;
					else singletonsV2++;
					line.str("");
					if(bipartiteOriginalEntities.size()>0)                                     
						line << "Community " << cont++ << "[" << _communities[i].getDescription() << "]: " << bipartiteOriginalEntities[g._graph[_communities[i].getNodes()[0]].getNodes()[0].getIdInput()] << "\n";
					else
						line << "Community " << cont++ << "[" << _communities[i].getDescription() << "]: " << g._graph[_communities[i].getNodes()[0]].getNodes()[0].getIdInput() << "\n";
						
					outfileC << line.str();
				}
				else
				{
					line.str("");
					line << "Community " << cont++ << "[" << _communities[i].getDescription() << "]: ";
					outfileC << line.str();
					line.str("");
					if(bipartiteOriginalEntities.size()>0)
					{
						for(int k=0;k<g._graph[_communities[i].getNodes()[0]].getNumberNodes();k++)
							line << bipartiteOriginalEntities[g._graph[_communities[i].getNodes()[0]].getNodesSorted()[k].getIdInput()] << ",";
					}
					else
					{
						for(int k=0;k<g._graph[_communities[i].getNodes()[0]].getNumberNodes();k++)
                                                        line << g._graph[_communities[i].getNodes()[0]].getNodesSorted()[k].getIdInput() << ",";
					}
					outfileC << line.str().substr(0,line.str().length()-1) << "\n";
				}
			}
			else
			{
				line.str("");
				line << "Community " << cont++ << "[" << _communities[i].getDescription() << "]: ";
				outfileC << line.str();
				line.str("");
				if(bipartiteOriginalEntities.size()>0)
				{
					for(int j=0;j<_communities[i].getNumberNodes();j++)
					{
						for(int k=0;k<g._graph[_communities[i].getNodes()[j]].getNumberNodes();k++)
							line << bipartiteOriginalEntities[g._graph[_communities[i].getNodes()[j]].getNodesSorted()[k].getIdInput()] << ",";
					}
				}
				else
				{
					for(int j=0;j<_communities[i].getNumberNodes();j++)
                                        {
                                                for(int k=0;k<g._graph[_communities[i].getNodes()[j]].getNumberNodes();k++)
                                                        line << g._graph[_communities[i].getNodes()[j]].getNodesSorted()[k].getIdInput() << ",";
                                        }

				}
				outfileC << line.str().substr(0,line.str().length()-1) << "\n";
			}
		}
	}
	line.str("");
	line << "\nSingletons Partition V1: " << singletonsV1;
	outfileC << line.str();
	line.str("");
	line << "\nSingletons Partition V2: " << singletonsV2 << "\n";
	outfileC << line.str();
	outfileC << listNodesCommunities(g);
	outfileC.close();
}


void biLouvainMethod::printCommunityNodes(int communityId)
{
	bool band = false;
	for(int i=0;i<_numberCommunities;i++)
	{
		if(_communities[i].getId()==communityId)
		{
			band = true;
			if (!_communities[i].getNodes().empty())
			{
				printf("\nCommunity ID: %d \n", _communities[i].getId());
				for (int j = 0; j < _communities[i].getNumberNodes(); j++)
					printf("\nNode: %d",_communities[i].getNodes()[j]);
			}
			break;
		}
	}
	if (band==false)
		printf("\nCommunity doesn't exist \n");
}

void biLouvainMethod::printCommunityNodesNeighbors(Graph &g,int communityId)
{
	bool band = false;
	for(int i=0;i<_numberCommunities;i++)
	{
		if(_communities[i].getId()==communityId)
		{
			band = true;
			if (~_communities[i].getNodes().empty())
			{
				for (int j = 0; j < _communities[i].getNumberNodes(); j++)
				{
					printf("\nCommunity ID: %d Node: %d \n", _communities[i].getId(),_communities[i].getNodes()[j]);
					if(g._graph[_communities[i].getNodes()[j]].getNeighbors().empty())
						printf("Node %d doesn't have neighbors \n", communityId);
					else
					{
						for (int k = 0; k <g. _graph[_communities[i].getNodes()[j]].getNumberNeighbors(); k++)
							printf("Neighbor: %d \n", g._graph[_communities[i].getNodes()[j]].getNeighbors()[k]);
					}
				}
			}
			break;
		}
	}
	if (band==false)
		printf("\nCommunity doesn't exist \n");
}

void biLouvainMethod::printTimes(double biLouvainTime, double loadGraphTime, double fusingTime)
{
	std::string outputfileTime = _outputFileName + "_ResultsTime.txt";
	std::ofstream outfileTime;
	outfileTime.open(outputfileTime.c_str(),std::ios::out|std::ios::trunc);
	outfileTime << "::: Total Time: " << timeConverter(biLouvainTime+loadGraphTime+fusingTime).c_str() << "microseconds: " << biLouvainTime+loadGraphTime+fusingTime << "\n";
	outfileTime << "\n::: Load Graph Total Time: " << timeConverter(loadGraphTime).c_str() << "microseconds: " << loadGraphTime << "\n";
	outfileTime << "\n::: Fuse Total Time: " << timeConverter(fusingTime).c_str() << "microseconds: " << fusingTime << "\n";
	outfileTime << "\n::: biLouvain Algorithm Total Time: " << timeConverter(biLouvainTime).c_str() << "microseconds: " << biLouvainTime << "\n";
	outfileTime << "\n::: Initial Times :::";
	outfileTime << "\nInitial Community Time: " + timeConverter(initialCommunityTime);
	outfileTime << "\nInitial Community Neighbors Time: " + timeConverter(initialCommunityNeighborsTime);
	outfileTime << "\nInitial CoCluster Time: " + timeConverter(initialCoClusterMateTime);
	outfileTime << "\n\n::: biLouvain Times :::";
	outfileTime << "\nbiLouvain Targets Time: " + timeConverter(candidatesTime);
	outfileTime << "\nbiLouvain Modularity Gain Time: " + timeConverter(modularityGainTime);
	outfileTime << "\nbiLouvain Update Time: " + timeConverter(updateTime);
	outfileTime << "\n\n::: biLouvain Gain PreCalculation Times :::";
	outfileTime << "\nbiLouvain Gain Time PC Ci: " + timeConverter(precalculationCiTime);
	outfileTime << "\nbiLouvain Gain Time PC Cj: " + timeConverter(precalculationCjTime);
	outfileTime << "\nbiLouvain Gain Time PC D: " + timeConverter(precalculationDTime);
	outfileTime << "\nPre Murata Time: " + timeConverter(premurataTime);
	outfileTime.close();
	/*printf("\n\n ::: Total Time: %s ::: %f microseconds\n",timeConverter(biLouvainTime+loadGraphTime).c_str(),biLouvainTime+loadGraphTime);
	printf("\n\n ::: Load Graph Total Time: %s ::: %f microseconds\n",timeConverter(loadGraphTime).c_str(),loadGraphTime);
	printf("\n\n ::: biLouvain Algorithm Total Time: %s ::: %f microseconds\n",timeConverter(loadGraphTime).c_str(),loadGraphTime);
	printf("\n::: Initial Times :::");
	printf("\nInitial Community Time: %s",timeConverter(initialCommunityTime).c_str());
	printf("\nInitial Community Neighbors Time: %s",timeConverter(initialCommunityNeighborsTime).c_str());
	printf("\nInitial CoCluster Time: %s",timeConverter(initialCoClusterMateTime).c_str());
	printf("\n\n::: biLouvain Times :::");
	printf("\nbiLouvain Targets Time: %s",timeConverter(candidatesTime).c_str());
	printf("\nbiLouvain Modularity Gain Time: %s",timeConverter(modularityGainTime).c_str());
	printf("\nbiLouvain Update Time: %s \n",timeConverter(updateTime).c_str());
	printf("\n::: biLouvain Gain PreCalculation Times :::");
	printf("\nbiLouvain Gain Time PC Ci: %s",timeConverter(precalculationCiTime).c_str());
	printf("\nbiLouvain Gain Time PC Cj: %s",timeConverter(precalculationCjTime).c_str());
	printf("\nbiLouvain Gain Time PC D: %s",timeConverter(precalculationDTime).c_str());
	printf("\nPre Murata Time: %s\n",timeConverter(premurataTime).c_str());*/
}
