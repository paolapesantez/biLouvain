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
# PreProcessInputBipartiteGraph.h
# Given the input file of the graph, we will generate a extra file with the required format to run biLouvain.
*/

#ifndef PREPROCESSINPUTBIPARTITEGRAPH_H_
#define PREPROCESSINPUTBIPARTITEGRAPH_H_

#include "Header.h"
#include "StringSplitter.h"

class PreProcessInputBipartiteGraph
{
	public:

	std::tr1::unordered_map<int,std::string> static preProcessingGraphData(const std::string &inputFileName, const std::string &delimiter)
	{
		std::tr1::unordered_map<int,std::string> bipartiteOriginalEntities;
		int pos = inputFileName.find_last_of(".");
                std::string bipartiteFileName = inputFileName.substr(0,pos)+"_bipartite.txt";
		std::string dictionaryFileName = inputFileName.substr(0,pos)+"_bipartite_Dictionary.txt";
                std::ifstream bipartiteFile(bipartiteFileName.c_str());
		if(bipartiteFile.is_open()==false)
                {
                	printf("\n ::: Reading and Preprocessing Input File...%s  Please Wait... :::",inputFileName.c_str());
			std::ifstream inputFile(inputFileName.c_str());
			std::string line = "";
			while(inputFile.good())
                        {
                        	getline(inputFile,line);
                              	int location = line.find(delimiter);
                        	if(location == std::string::npos) 
				{
					printf("\n ::: Given delimiter doesn't agree with delimiter being used in file :::\n");
	                                exit(EXIT_FAILURE);
				}
				break;
			}
			inputFile.clear();
                        inputFile.seekg(0, std::ios::beg);					
			std::tr1::unordered_map<std::string,int> bipartiteIds;
			std::ofstream bipartiteOFile;
			std::ofstream dictionaryFile;
			dictionaryFile.open(dictionaryFileName.c_str(),std::ios::out|std::ios::trunc);
			bipartiteOFile.open(bipartiteFileName.c_str(),std::ios::out|std::ios::trunc);
			int cont = 0;
			std::string *pieces;
                        int items = 0;
			std::stringstream entry;	
			while(inputFile.good())
			{
				getline(inputFile,line);
				if(inputFile.eof())break;
				pieces = StringSplitter::split(line,delimiter,items);
				//std::cout<<"\n" << pieces[0] <<"   " << items;
				if(bipartiteIds.find(pieces[0].c_str()) == bipartiteIds.end())
				//if(bipartiteOriginalEntities.end() == find_if(bipartiteOriginalEntities.begin(),bipartiteOriginalEntities.end(),[&value](const map_value_type& vt){ return vt.second == pieces[0]}))
				{
					bipartiteIds[pieces[0]] = cont;	
					bipartiteOriginalEntities[cont] = pieces[0];
					entry.str("");
					entry << cont << "\t" << pieces[0] << "\n";
					dictionaryFile << entry.str();
					cont++;
				}
			}
			inputFile.clear();
			inputFile.seekg(0,std::ios::beg);
			while(inputFile.good())
			{
				std::stringstream edge;
				getline(inputFile,line);
				if(inputFile.eof())break;
				pieces = StringSplitter::split(line,delimiter,items);
				//std::cout<<"\n" << pieces[1];
				//if(bipartiteOriginalEntities.end() == find_if(bipartiteOriginalEntities.begin(),bipartiteOriginalEntities.end(),[&value](const map_value_type& vt){ return vt.second == pieces[1]}))
				if(bipartiteIds.find(pieces[1].c_str()) == bipartiteIds.end())
				{
					bipartiteIds[pieces[1]] = cont;
					bipartiteOriginalEntities[cont] = pieces[1];
					entry.str("");
					entry << cont << "\t" << pieces[1] << "\n";
					dictionaryFile << entry.str();
					cont++;
				}
				edge.str("");
				//edge.precision(15);
				if(items == 3)
					edge << bipartiteIds[pieces[0]] << "\t" << bipartiteIds[pieces[1]] << "\t" << pieces[2] << "\n";
				else
					edge << bipartiteIds[pieces[0]] << "\t" << bipartiteIds[pieces[1]] << "\t1" << "\n";
				bipartiteOFile << edge.str();
			}
			inputFile.close();
			dictionaryFile.close();
			bipartiteOFile.close();
			bipartiteIds.clear();
			delete[] pieces;	
		}
		else
			bipartiteOriginalEntities = readDictionaryFile(bipartiteFileName);
		return bipartiteOriginalEntities;
	}

	
	std::tr1::unordered_map<int,std::string> static readDictionaryFile(const std::string &inputFileName)
       	{
                std::tr1::unordered_map<int,std::string> bipartiteOriginalEntities;
                int pos = inputFileName.find_last_of(".");
                std::string dictionaryFileName = inputFileName.substr(0,pos)+"_Dictionary.txt";
		std::ifstream bipartiteDictionaryFile(dictionaryFileName.c_str());
		if(bipartiteDictionaryFile.is_open()==true)
		{
			printf("\n ::: Reading Dictionary File...%s  Please Wait... :::",inputFileName.c_str());
                        std::string line = "";
	        	std::string *pieces;
        	        int items = 0;
                        while(bipartiteDictionaryFile.good())
	                {
        			getline(bipartiteDictionaryFile,line);
                	        if(bipartiteDictionaryFile.eof())break;
                                pieces = StringSplitter::split(line,"\t",items);
	                        //std::cout<<"\n" << pieces[0] <<"   " << items;
                                bipartiteOriginalEntities[stoi(pieces[0])] = pieces[1];
                        }
			bipartiteDictionaryFile.close();
                        delete[] pieces;
                }
		else
			printf("\n ::: Warning: Dictionary File was not found. :::");			
		return bipartiteOriginalEntities;
	}
};
#endif
