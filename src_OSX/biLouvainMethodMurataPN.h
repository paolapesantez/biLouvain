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
# biLouvainMethodMurataPN.h
# It contains the main methods to calculate Murata+ Modularity
*/

#ifndef BILOUVAINMETHODMURATAPN_H_
#define BILOUVAINMETHODMURATAPN_H_


#include "biLouvainMethod.h"

class biLouvainMethodMurataPN : public biLouvainMethod
{
	protected:
		int calculateNumberNodesBetaFactor(Graph &g,int &communityId, MetaNode &node, int option);
		double calculateCommunityBetaFactor(Graph &g,std::string communityType,double similarity);
		double calculateCommunitySimilarity(Graph &g,int &communityId);
		double murataModularityArgMax(Graph &g,int &communityId, int possibleCoClusterMateId);
		double* murataModularityWithChanges(Graph &g,MetaNode &node, int &communityId, int possibleCoClusterMateId,int &newCommunityId,int &option);
		newDataCommunityVector murataCalculationCoClusterMates(Graph &g,int communityId,const std::vector<int> &possibleCoClusterMates);
		double CoClusterMateDefinitionAllCommunities(Graph &g, int start, int end);
		double IntraTypeDefinitionAllCommunities(Graph &g, int start, int end);
		void CoClusterMateDefinitionIDCommunity(Graph &g,int &communityId);
		newDataCommunity CoClusterMateDefinitionPrecalculation(Graph &g,MetaNode &node, int &communityId, int &newCommunityId,int &optio);

	public:
		biLouvainMethodMurataPN();
		~biLouvainMethodMurataPN();

};

#endif /* BILOUVAINMETHODMURATAPN_H_ */
