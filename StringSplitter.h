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
# StringSplitter.h
# Parses a line according a delimiter, dividing the line in its different parts. 
# It also provides methods to print the element of any kind of vector or array.
*/

#ifndef STRINGSPLITTER_H_
#define STRINGSPLITTER_H_


#include "Header.h"


class StringSplitter
{
	public:

		//Accepts a std::string and a delimiter.  Will use items_found to return the number
		//of items found as well as an array of std::strings where each element is a piece of
		//the original std::string.
		static std::string * split(std::string text, std::string delimiter, int &items_found)
		{
			//std::vectors are dynamically expanding arrays
			std::vector<std::string> pieces;

			//find the first delimiter
			int location = text.find(delimiter);

			//we are starting at the beginning of our std::string
			int start = 0;

			//go until we have no more delimiters
			while(location != std::string::npos)
			{
				//add the current piece to our list of pieces
				std::string piece = text.substr(start, location - start);
				pieces.push_back(piece);

				//update our index markers for the next round
				start = location + 1;
				location = text.find(delimiter, start);
			}

			//at the end of our loop, we're going to have one trailing piece to take care of.
			//handle that now.
			std::string piece = text.substr(start, location - start);
			pieces.push_back(piece);

			//convert from std::vector into an array of std::strings
			int size = pieces.size();
			std::string *pieces_str = new std::string[size];
			for(int i = 0; i < size; i++)
			{
				pieces_str[i] = pieces.at(i);
			}
			items_found = size;
			return pieces_str;
		}

		template <typename T>
		static void printVector(std::vector<T> _vector)
		{
			if(_vector.size()>0)
			{
				std::stringstream element;
				for(unsigned int i=0;i<_vector.size();i++)
				{
					element << _vector[i];
					printf("\nElement: %s", (element.str()).c_str());
					element.str("");
				}
			}
			else
				std::cout << "Vector is empty!" << std::endl;
		}

		template <typename T>
		static void printArray(T*_array,int n)
		{
			std::stringstream element;
			for(int i=0;i<n;i++)
			{
				element << _array[i];
				printf("\nElement: %s", (element.str()).c_str());
				element.str("");
			}
		}
};

#endif /* STRINGSPLITTER_H_ */
