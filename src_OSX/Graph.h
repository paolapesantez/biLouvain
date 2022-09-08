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
# Graph.h
# 
*/

#ifndef GRAPH_H_
#define GRAPH_H_

#include "Header.h"
#include "MetaNode.h"

class Graph
{
	friend class biLouvainMethod;
	friend class biLouvainMethodMurataPN;
	friend class FuseMethod;

	protected:
		MetaNode* _graph;
		int _numberNodes;
		int _numberEdges;
		double _weightEdges;
		double _weightEdgesV1;
		double _weightEdgesV2;
		int _lastIdPartitionV1;
		
		double _lambdaV1;
		double _lambdaV2;

		double _sumSimilarityV1;
		double _sumSimilarityV2;

	public:
		//Class constructor
		Graph();
		Graph(MetaNode* &graph,int &numberNodes,int &numberEdges,double &weightEdges,double &weightEdgesV1,double &weightEdgesV2,int &lastIdPartitionV1);
		//Class destructor
		~Graph();
		MetaNode* getGraph();
		MetaNode getNode(const int &nodeId);
		int getNumberNodes();
		int getLastIdPartitionV1();
		int getNumberEdges();
		double getWeightEdges();
		int getNumberEdgesV1();
		int getNumberEdgesV2();
		double getWeightEdgesV1();
		double getWeightEdgesV2();
		double getLambdaV1();
		double getLambdaV2();
		double getSimilarityV1();
		double getSimilarityV2();
		void setLambdaV1(double lambdaV1);
		void setLambdaV2(double lambdaV2);
		void setSimilarityV1(double sumSimilarityV1);
		void setSimilarityV2(double sumSImilarityV2);
		void setNumberCoClusters(int number);
		void addIntraTypeNeighborsToNode(int &nodeId,std::unordered_map<int,double> &intraTypeNeighbors);
		void printNeighborsNode(int nodeId);
		void printGraph(const std::string &inputFileName);
		void destroyGraph();

};

#endif /* GRAPH_H_ */
