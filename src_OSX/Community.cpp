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


#include "Community.h"

Community::Community(){}

Community::Community(int id, std::string description, std::unordered_map<int,double> nodes,int option)
{
        _id   = id;
        _description = description;
	if(option == 1)
	        _nodes = nodes;
	else
		_nodesIntraType = nodes;
        _modularityContribution = 0.0;
        _betaFactor = 0.0;
}


Community::Community(int id, std::string description, std::unordered_map<int,double> nodes,std::unordered_map<int,double> nodesIntraType)
{
	_id   = id;
	_description = description;
	_nodes = nodes;
	_nodesIntraType = nodesIntraType;
	_modularityContribution = 0.0;
	_betaFactor = 0.0;
}

/* Get functions */
int Community::getId()
{
	return _id;
}

std::string Community::getDescription()
{
	return _description;
}

double Community::getModularityContribution()
{
	return _modularityContribution;
}

double Community::getBetaFactor()
{
        return _betaFactor;
}

std::vector<int> Community::getNodes()
{
	std::vector<int> result;
	if(_nodes.size() != 0)
	{
		for(auto it=_nodes.begin();it!=_nodes.end();++it)
			result.push_back(it->first);
	}
	else
	{
		for(auto it=_nodesIntraType.begin();it!=_nodesIntraType.end();++it)
                        result.push_back(it->first);	
	}
	return result;
}

std::vector<int> Community::getCoClusterMateCommunityId()
{
	std::vector<int> result;
	result = _coClusterMateCommunityId;
	return result;
}

int Community::getNumberNodes()
{
	if(_nodes.size() != 0)
		return _nodes.size();
	else
		return _nodesIntraType.size();
}

double Community::getDegreeCommunity()
{
	double result = 0.0;
	for(auto it=_nodes.begin();it!=_nodes.end();++it)
		result += it->second;
	return result;
}

double Community::getDegreeCommunityWithoutNode(int node_id)
{
	double result = 0.0;
	for(auto it=_nodes.begin();it!=_nodes.end();++it)
		if(it->first != node_id)
			result += it->second;
	return result;
}

double Community::getSimilarity()
{
	double result = 0.0;
        for(auto it=_nodesIntraType.begin();it!=_nodesIntraType.end();++it)
        	result += it->second;
        return result;
}

double Community::getSimilarityWithoutNode(int nodeId)
{
        double result = 0.0;
        for(auto it=_nodesIntraType.begin();it!=_nodesIntraType.end();++it)
                if(it->first != nodeId)
                        result += it->second;
        return result;
}

std::vector<int> Community::getNodesWithoutNode(int nodeId)
{
	std::vector<int> result;
	for(auto it=_nodes.begin();it!=_nodes.end();++it)
		if(it->first != nodeId) 
			result.push_back(it->first);
	return result;
}

/* Set procedures */
void Community::setId(int id)
{
	_id = id;
}

void Community::setDescription(std::string description)
{
	_description = description;
}

void Community::setCoClusterMateCommunityId(std::vector<int> coClusterMateCommunityId)
{
	_coClusterMateCommunityId.clear();
	_coClusterMateCommunityId = coClusterMateCommunityId;
}

void Community::setNodes(std::unordered_map<int,double> nodes)
{
	_nodes = nodes;
}

void Community::setNodesIntratype(std::unordered_map<int,double> nodesIntraType)
{
        _nodesIntraType = nodesIntraType;
}

void Community::setModularityContribution(double modularityContribution)
{
	_modularityContribution = modularityContribution;
}

void Community::setBetaFactor(double betaFactor)
{
        _betaFactor = betaFactor;
}

void Community::addNode(int nodeId, double nodeDegree)
{
	_nodes[nodeId] = nodeDegree;
}


void Community::deleteNode(int nodeId)
{
	_nodes.erase(nodeId);
}

void Community::addIntraTypeNode(int nodeId, double nodeSimilarity)
{
        _nodesIntraType[nodeId] = nodeSimilarity;
}


void Community::deleteIntraTypeNode(int nodeId)
{
        _nodesIntraType.erase(nodeId);
}
