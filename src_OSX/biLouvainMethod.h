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


/*
# biLouvainMethod.h
# It contains the main methods for the execution of the biLouvain Algorithm.
*/

#ifndef BILOUVAINMETHOD_H_
#define BILOUVAINMETHOD_H_

#include "StringSplitter.h"
#include "Graph.h"
#include "MetaNode.h"
#include "Node.h"
#include "Community.h"
#include "Timer.h"

class biLouvainMethod
{
	protected:
		double _alpha;	
		double _totalPartitioningModularity;	
		int _numberCommunities;
		int _numberCommunitesV1;
		int _numberCommunitiesV2;
		std::vector<Community> _communities;
		std::string _outputFileName;
		//std::vector<double> _communitiesBetaFactor;
		std::vector<Node> nodes;

		/*Auxiliar Functions and Procedures*/
		int findCommunityContainingNode(int nodeId);
		int isRepeated(const std::vector<int> &elements, int newElement);
		int calculateEdgesBetweenCommunities(Graph &g,int communityCId, int communityDId);
		int findPositionNode(Graph &g,int node_id);
		double calculateEdgesBetweenCommunitiesMap(Graph &g,int communityCId, int communityDId);
		std::vector<int> findNeighborCommunities(Graph &g,int communityId);
		std::vector<int> findNeighborCommunitiesMap(Graph &g,int communityId);
		std::vector<int> findNeighborCommunitiesWithoutNodeMap(Graph &g,int communityId, int nodeId);
		std::vector<int> getDifferentNeighborCommunities(Graph &g,int communityId1, int communityId2);
		std::vector<int> getDifferentNeighborCommunitiesMap(Graph &g,int communityId1, int communityId2);
		void calculateCommunitiesModulatiryContribution();

		/*Main Functions and Procedures*/
		void initialCommunityDefinition(Graph &g);
		void initialCommunityDefinitionIntraType(Graph &g);
		void initialCommunityDefinitionWithIntraType(Graph &g);
		void initialCommunityNeighborsDefinition(Graph &g);
		void initialIntraTypeCommunityNeighborsDefinition(Graph &g);
		void updateNodeCommunity(Graph &g,int nodeId, int oldCommunityId, int newCommunityId);
		void updateNodeIntraTypeCommunity(Graph &g,int nodeId, int oldCommunityId, int newCommunityId);
		void updateCoClusterMateCommunities(const std::string &changes);
		void updateNeighborCommunities(Graph &g,int nodeId, int oldCommunityId,int newCommunityId);
		void updateIntraTypeNeighborCommunities(Graph &g,int nodeId, int oldCommunityId,int newCommunityId);
		std::unordered_map<int,double> compactMetaNodeNeighbors(Graph &g,int &communityId, std::unordered_map<int,int> &dictionaryCommunities);
		std::unordered_map<int,double> compactMetaNodeIntraTypeNeighbors(Graph &g,int &communityId, std::unordered_map<int,int> &dictionaryCommunities);
		void fromCommunitiesToNodes(Graph &g);
		int numberCommunitiesNonEmpty();
		std::unordered_map<int,int> dictionaryCommunitiesNewId();
		virtual double murataModularityArgMax(Graph &g,int &communityId, int possibleCoClusterMateId)=0;
		virtual double* murataModularityWithChanges(Graph &g,MetaNode &node, int &communityId, int possibleCoClusterMateId,int &newCommunityId,int &option)=0;
		virtual double CoClusterMateDefinitionAllCommunities(Graph &g,int start, int end)=0;
		virtual double IntraTypeDefinitionAllCommunities(Graph &g,int start, int end)=0;
		virtual newDataCommunity CoClusterMateDefinitionPrecalculation(Graph &g,MetaNode &node, int &communityId, int &newCommunityId,int &option)=0;
		virtual double calculateCommunityBetaFactor(Graph &g,std::string communityType,double similarity)=0;
                virtual double calculateCommunitySimilarity(Graph &g,int &communityId)=0;
		double calculateMaxModularityGainIteration(Graph &g,int* &nodesOrderExecution);
		double calculateMaxModularityGainIterationIntraType(Graph &g,int* &nodesOrderExecution);
		newDataCommunity calculateDeltaGainModularity(Graph &g,MetaNode &node, int &communityId, int newCommunityId,int option);
		int* nodesOrderToProcess(Graph &g,int optionOrder);

	public:
		double initialCommunityTime;
		double initialCommunityNeighborsTime;
		double initialCoClusterMateTime;
		double candidatesTime;
		double modularityGainTime;
		double updateTime;
		double precalculationCiTime;
		double precalculationCjTime;
		double precalculationDTime;
		double premurataTime;

		biLouvainMethod();
		~biLouvainMethod();
		/*biLouvain Method*/
		void biLouvainMethodAlgorithm(Graph &g,double cutoffIterations, double cutoffPhase, int optionOrder,std::unordered_map<int,std::string> &bipartiteOriginalEntities,const std::string &inputFileName,const std::string &outputFileName,double &alpha);
		void biLouvainMethodAlgorithmIntraType(Graph &g,double cutoffIterations, double cutoffPhase, int optionOrder,std::unordered_map<int,std::string> &bipartiteOriginalEntities,const std::string &inputFileName,const std::string &outputFileName);

		/*Printing and storing results*/
		int numberNodesInsideCommunity(Graph &g,int communityId);
		std::string listNodesCommunities(Graph &g);
		void generateOutputFile(std::string text, std::string fileName);
		void printCommunitiesContributionModularity();
		void printCommunities(Graph &g);
		void printCoClusterCommunitiesFile();
		void printAllCommunityNodes(Graph &g);
		void printAllCommunityNodeswithSingletons(Graph &g,std::unordered_map<int,std::string> &bipartiteOriginalEntities);
		void printCommunityNodes(int communityId);
		void printCommunityNodesNeighbors(Graph &g,int communityId);
		void printTimes(double totalTime, double loadGraphTime,double mergingTime);
};

#endif /* BILOUVAINMETHOD_H_ */
