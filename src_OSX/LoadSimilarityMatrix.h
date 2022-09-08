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
# LoadSimilarityMatrix.h
# Reads the input file/files containing the similarity matrices of each set of vertices of the bipartite graph.
*/

#ifndef LOADSIMILARITYMATRIX_H_
#define LOADSIMILARITYMATRIX_H_
#include "Header.h"
#include "StringSplitter.h"

class LoadSimilarityMatrix
{
  public:


	//Read the similarity matrices for V1, V2, or both.
	int static loadSimilarityMatrixFromFile(Graph &g,std::string &similarityMatrixFileName)
	{
		std::ifstream inputFile(similarityMatrixFileName.c_str());
		bool is_good = true;
		int result;

		if(inputFile.is_open() == false)			//If the file hasn't been found
		{
			result = -1;							//no file has been found. End the program.
			is_good = false;
		}
		if(is_good)						//If the file has been found
		{
			std::string line = "";
			getline(inputFile,line);
			while(line[0] == '#')
				getline(inputFile,line);
	                int location = line.find(":");
 			if(location == std::string::npos)
                        {
                             	printf("\n ::: Expected similarity (alpha) not specified :::\n");
                             	exit(EXIT_FAILURE);
                        }
                        else
			{
				int row = 0, column = 0;
				double lambdaPartition =0.0, sumSimilarityPartition = 0.0;
				lambdaPartition = atof(line.substr(location+1,line.length()).c_str());
				if(similarityMatrixFileName.find("V2") != std::string::npos)
                                        column = g.getLastIdPartitionV1()+1;
				getline(inputFile,line);
				location = line.find(",");
				if(location == std::string::npos)
				{
					printf("\n ::: The similarity matrix should be "," delimited :::\n");
					exit(EXIT_FAILURE);
				}
				inputFile.clear();
	                        inputFile.seekg(0, std::ios::beg);
				int columnId = 0, rowId = 0;
				int items = 0;
				std::string* pieces;
				std::unordered_map<int,double> intraTypeNeighbors;
				while(inputFile.good())				//read line by line
				{
					getline(inputFile,line);
					//std::cout << "\n" << line;
					if(inputFile.eof())break;
					if((line.length()>0)&&(line[0] != '#')&&(line.find(":") == std::string::npos))
					{       
						pieces = StringSplitter::split(line,",",items);
						for(int i=0; i<items;i++)						
						{
							rowId = row + column;
							columnId = column + i;
							if((atof(pieces[i].c_str()) != 0.0)&&(rowId != columnId))
							{
							    intraTypeNeighbors[columnId] = atof(pieces[i].c_str());
 							    sumSimilarityPartition += atof(pieces[i].c_str());
							}
						}
						g.addIntraTypeNeighborsToNode(rowId,intraTypeNeighbors);
						intraTypeNeighbors.clear();
						row++;
					}
				}		
				std::cout<<"\n Similarity Sum:"<< sumSimilarityPartition;
				if(similarityMatrixFileName.find("V1") != std::string::npos)
				{
                                        g.setSimilarityV1(sumSimilarityPartition);
					g.setLambdaV1(lambdaPartition/sumSimilarityPartition);
				}
                                else	
				{
					g.setSimilarityV2(sumSimilarityPartition);
					g.setLambdaV2(lambdaPartition/sumSimilarityPartition);
				}
				delete[] pieces;
				inputFile.close();		//close the file from which we were reading
				result = 0;				//reading completed successfully
			}					
		}
		return result;
	}

};

#endif /* LOADSIMILARITYMATRIX_H_ */
