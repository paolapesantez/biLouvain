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
# LoadGraph.h
# Reads the input file indicated by the user and preprocesses it into the format needed for biLouvain to work.
*/

#ifndef LOADGRAPH_H_
#define LOADGRAPH_H_
#include "StringSplitter.h"
#include "Node.h"
#include "MetaNode.h"

class LoadGraph
{
  public:


	//Read the bipartite graph from the input file input by the user
	int static loadBipartiteGraphFromFile(Graph* &g,std::string &inputFileName)
	{
		MetaNode* _graph = NULL;
		int _numberNodes = 0;
		int _numberEdges = 0;
		double _weightEdges = 0.0;
		double _weightEdgesV1 = 0.0;
		double _weightEdgesV2 = 0.0;
		int _lastIdPartitionV1 = 0;
		std::ifstream inputFile(inputFileName.c_str());
		bool is_good = true;
		int result;

		if(inputFile.is_open() == false)			//If the file hasn't been found
		{
			result = -1;					//no file has been found. End the program.
			is_good = false;
		}
		if(is_good)						//If the file has been found
		{
			std::string line = "";
                        getline(inputFile,line);
                        int location = line.find("\t");
                        if(location == std::string::npos)
                        {
                                printf("\n ::: The bipartite file should be tab delimited :::\n");
                        	exit(EXIT_FAILURE);
			}
			//Getting the maximum id which will indicate us the number of nodes in the bipartite graph
			inputFile.clear();
                        inputFile.seekg(0, std::ios::beg);
			int maximumId = 0;
			int items = 0;
			std::string* pieces;
			while(inputFile.good())				//read line by line
			{
				getline(inputFile,line);
				if(inputFile.eof())break;
				if((line.length()>0)&&(line[0] != '#')) //because some inputs can contain informative lines starting with #
				{
					pieces = StringSplitter::split(line,"\t",items);
					if(atoi(pieces[1].c_str())>maximumId)
						maximumId = atoi(pieces[1].c_str());
					if(atoi(pieces[0].c_str())>_lastIdPartitionV1)
						_lastIdPartitionV1 = atoi(pieces[0].c_str());
				}
			}
			_numberNodes = maximumId+1;
			_graph = new MetaNode[_numberNodes];
			inputFile.clear();
			inputFile.seekg(0, std::ios::beg);
	
			std::unordered_map<int,double> neighborsPerNode[_numberNodes];
			while(inputFile.good())				//read line by line
			{
				getline(inputFile,line);
				if(inputFile.eof())break;
				if((line.length()>0)&&(line[0] != '#'))
				{       //find neighbors for each node
					pieces = StringSplitter::split(line,"\t",items);
					neighborsPerNode[atoi(pieces[0].c_str())][atoi(pieces[1].c_str())] = atof(pieces[2].c_str());		
					neighborsPerNode[atoi(pieces[1].c_str())][atoi(pieces[0].c_str())] = atof(pieces[2].c_str());
					_numberEdges++;
					_weightEdges += atof(pieces[2].c_str());
				}
			}
			std::vector<Node> nodeV;
			for(int i=0;i<_numberNodes;i++)   //Creation of metanodes and nodes
			{					
				if(i<=_lastIdPartitionV1)	//Create nodes belonging to set V1
				{	Node node(i,"V1",0);
					nodeV.push_back(node);
					MetaNode metanode(i,"V1",nodeV,neighborsPerNode[i],-1);
					for(auto it=neighborsPerNode[i].begin();it!=neighborsPerNode[i].end();++it)
						_weightEdgesV1 += it->second;
					_graph[i] = metanode;
					nodeV.clear();
				}
				else				//Create nodes belonging to set V2
				{	Node node(i,"V2",0);	
					nodeV.push_back(node);
					MetaNode metanode(i,"V2",nodeV,neighborsPerNode[i],-1);
					for(auto it=neighborsPerNode[i].begin();it!=neighborsPerNode[i].end();++it)
						_weightEdgesV2 += it->second;
					_graph[i] = metanode;
					nodeV.clear();
				}
			}
			delete[] pieces;
			inputFile.close();		//close the file from which we were reading
			result = 0;			//reading completed successfully
			g = new Graph(_graph,_numberNodes,_numberEdges,_weightEdges,_weightEdgesV1,_weightEdgesV2,_lastIdPartitionV1);
			//delete[] _graph;
		}
		return result;
	}

};

#endif /* LOADGRAPH_H_ */
