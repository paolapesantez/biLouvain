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

#include "Graph.h"

//Class constructor
Graph::Graph()
		:_graph(NULL),_numberNodes(0),_numberEdges(0),_weightEdges(0.0),_weightEdgesV1(0.0),_weightEdgesV2(0.0),_lastIdPartitionV1(0)
{}

Graph::Graph(MetaNode* &graph,int &numberNodes,int &numberEdges,double &weightEdges,double &weightEdgesV1,double &weightEdgesV2,int &lastIdPartitionV1)
		:_graph(graph),_numberNodes(numberNodes),_numberEdges(numberEdges),_weightEdges(weightEdges),_weightEdgesV1(weightEdgesV1),_weightEdgesV2(weightEdgesV2),_lastIdPartitionV1(lastIdPartitionV1)
{}

//Class destructor
Graph::~Graph(){}

//Read a file that contains info of a unipartite graph and turn it into a bipartite graph

MetaNode* Graph::getGraph()
{
	MetaNode* result;
	result = _graph;
	return result;
}

MetaNode Graph::getNode(const int &nodeId)
{
	return _graph[nodeId];		
}

int Graph::getNumberNodes()
{
	return _numberNodes;
}

int Graph::getLastIdPartitionV1()
{
	return _lastIdPartitionV1;
}

int Graph::getNumberEdges()
{
	return _numberEdges;
}

double Graph::getWeightEdges()
{
	return _weightEdges;
}


double Graph::getWeightEdgesV1()
{
	return _weightEdgesV1;
}

double Graph::getWeightEdgesV2()
{
	return _weightEdgesV2;
}


void Graph::printNeighborsNode(int id)
{
	for(int i=0;i<_numberNodes;i++)
	{
		if(_graph[i].getId()==id)
		{
			if(_graph[i].getNeighbors().empty())
			{
				printf("Node %d doesn't have neighbors \n", id);
			}
			else
			{
				printf("\nNode ID: %d \n", _graph[i].getId());
				for (int j = 0; j < _graph[i].getNumberNeighbors(); j++)
				{
					printf("Neighbor ID: %d \n",_graph[i].getNeighbors()[j]);
				}
			}
			break;
		}
	}
}

void Graph::printGraph(const std::string &inputFileName)
{
        std::string outputGraph = "FinalGraph" + inputFileName;
        std::ofstream outfileGraph;
        outfileGraph.open(outputGraph.c_str(),std::ios::out|std::ios::trunc);
        std::stringstream line;
        int i=0;
        while(_graph[i].getType()=="V1")
        {
                for(int j=0;j<_graph[i].getNumberNeighbors();j++)
                {
                        line.str("");
                        int neighbor = _graph[i].getNeighborsSorted()[j];
                        line << _graph[i].getId() << "\t" << neighbor << "\t" <<  _graph[i].getWeightNeighbor(neighbor) << std::endl;
                        outfileGraph << line.str();
                }
                i++;
        }
        outfileGraph.close();
}


void Graph::destroyGraph()
{
	delete[] _graph;
}
