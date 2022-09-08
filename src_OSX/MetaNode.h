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
# MetaNode.h
# Class that represents a node after the compaction process.
*/


#ifndef METANODE_H_
#define METANODE_H_

#include "Header.h"
#include "Node.h"

class MetaNode

{
	private:
		int _idGraph;
		std::string _type;
		int _communityId;
		std::vector<Node>_nodes;
		std::unordered_map<int,double> _neighbors;
		std::unordered_map<int,double> _neighborCommunities;
		std::unordered_map<int,double> _intraTypeNeighbors;
		std::unordered_map<int,double> _intraTypeNeighborCommunities;	
		struct CompareById {
			bool operator()(Node i, Node j) {return (i.getIdInput() < j.getIdInput());}
		}myobject;
		struct CompareInt {
			bool operator()(int i, int j) {return (i < j);}
		}myobject2;


	public:
		MetaNode();
		MetaNode(int id, std::string type, std::vector<Node> nodes,std::unordered_map<int,double> neighbors,int communityId);

		/* Get functions */
		int getId();
		std::string getType();
		int getCommunityId();
		int getNumberNodes();
		std::vector<Node> getNodes();
		std::vector<Node> getNodesSorted();
		std::vector<int> getNeighbors();
		std::vector<int> getNeighborsSorted();
		std::vector<int> getNeighborsWithoutNode(int nodeId);
		std::vector<int> getNeighborCommunities();
                std::vector<int> getIntraTypeNeighbors();
                std::vector<int> getIntraTypeNeighborCommunities();

		int getNumberNeighbors();
		int getNumberNeighborCommunities();
		int getNumberIntraTypeNeighbors();
		int getNumberIntraTypeNeighborCommunities();
		
		double getDegreeNode();
		double getDegreeNeighborCommunities();
		double getWeightEdgesToNeighborCommunity(int communityId);
		double getWeightNeighbor(int id);
		double getSimilarityNode();
		double getSimilarityIntraTypeNeighbor(int id);
		double getSimilarityToIntraTypeNeighborCommunity(int communityId);


		/* Set procedures */
		void setId(int id);
		void setType(std::string type);
		void setNodes(std::vector<Node> nodes);
		void setNeighbors(std::unordered_map<int,double> neighbors);
		void setCommunityId(int communityId);
		void setNeighborCommunities(std::unordered_map<int,double> neighborCommunities);
		void setIntraTypeNeighbors(std::unordered_map<int,double> intraTypeNeighbors);
		void setIntraTypeNeighborCommunities(std::unordered_map<int,double> intraTypeNeighborCommunities);
		void deleteNeighborCommunity(int communityId);
		void deleteNeighborCommunityWeight(int communityId, double weight);
		void deleteIntraTypeNeighborCommunitySimilarity(int communityId, double similarity);
		void addNeighbor(int idNeighbor,double weight);
		void addNeighborCommunity(int communityId);
		void addIntraTypeNeighborCommunitySimilarity(int communityId, double similarity);
		void addNeighborCommunityWeight(int communityId, double weight);
		int findNeighborCommunity(int communityId);
};

#endif /* METANODE_H_ */
