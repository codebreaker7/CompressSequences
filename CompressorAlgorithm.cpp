// CompressorAlgorithm.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include "BitIoStream.hpp"
#include "CanonicalCode.hpp"
#include "FrequencyTable.hpp"
#include "HuffmanCoder.hpp"
#include "rlehuffmancompr.h"

long get_file_size(const char* filename) {
	struct stat statv;
	int res = stat(filename, &statv);
	return res == 0 ? statv.st_size : -1;
}

int main(int argc, char * argv[])
{
    std::ifstream infile;
	if (argc < 2)
	{
		std::cout << "Need to provide a file name" << std::endl;
		return -1;
	}
	// prepare data
	infile.open(argv[1], std::ifstream::in);
	std::vector<std::pair<int, int>> compressVector[3];
	int prevValue[] = { -1, -1, -1 };
	int repeatCounter[] = { 0, 0, 0 };
	int length128[] = { 0, 0, 0 };
	int lengthnon128[] = { 0, 0, 0 };
	std::pair<int, int> tempPair;
	unsigned int numLines = 0;
	while (!infile.eof())
	{
		numLines++;
		int temp = 0;
		for (size_t j = 0; j < 3; j++)
		{
			//infile >> std::hex >> temp;
			infile.read((char*)&temp, sizeof(unsigned char));
			if (temp == prevValue[j])
			{
				repeatCounter[j]++;
				if (repeatCounter[j] == 65536)
				{
					tempPair.first = prevValue[j];
					tempPair.second = 65535;
					compressVector[j].push_back(tempPair);
					if (temp == 128)
					{
						length128[j]++;
					}
					repeatCounter[j] = 1;
				}
			}
			else
			{
				if (repeatCounter[j] >= 1)
				{
					tempPair.first = prevValue[j];
					tempPair.second = repeatCounter[j];
					compressVector[j].push_back(tempPair);
					if (prevValue[j] == 128 && repeatCounter[j] > 4)
					{
						length128[j]++;
					}
				}
				prevValue[j] = temp;
				repeatCounter[j] = 1;
			}
		}
	}
	numLines--;
	// append the last chunks of data
	for (size_t j = 0; j < 3; j++)
	{
		tempPair.first = prevValue[j];
		tempPair.second = repeatCounter[j];
		compressVector[j].push_back(tempPair);
		if (prevValue[j] == 128 && repeatCounter[j] > 4)
		{
			length128[j]++;
		}
	}
	infile.close();

	FrequencyTable frequencyTable(std::vector<uint32_t>(256, 0));
	unsigned int offset = 0;
	if (argc == 3) // the last parameter is a name of distribution file
	{
		frequencyTable = loadFrequencyTable(argv[argc - 1]);
	} 
	else { // otherwise calculate our table
		// calculate total length of non-compressed sequences in each column
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < compressVector[i].size(); j++)
			{
				if (compressVector[i].at(j).first != 128)
				{
					lengthnon128[i] += compressVector[i].at(j).second;
					offset += compressVector[i].at(j).second;
					frequencyTable.set(compressVector[i].at(j).first, frequencyTable.get(compressVector[i].at(j).first) + compressVector[i].at(j).second);
				}
				else if (compressVector[i].at(j).first == 128 && compressVector[i].at(j).second <= 4) {
					lengthnon128[i] += compressVector[i].at(j).second;
					offset += compressVector[i].at(j).second;
					frequencyTable.set(compressVector[i].at(j).first, frequencyTable.get(compressVector[i].at(j).first) + compressVector[i].at(j).second);
				}
				else
				{
					if (offset > 65535)
					{
						length128[i] += ((offset / 65535) - 1);
						length128[i] += (offset % 65535 != 0) ? 1 : 0;
					}
					offset = 0;
				}
			}
		}
	}
	for (size_t i = 0; i < 256; i++)
	{
		std::cout << i << "-" << frequencyTable.get(i) << std::endl;
	}
	// store 128-based sequences
	std::ofstream text_compress("compress.txt", std::ios::out);
	std::ofstream file128;
	std::ofstream file_residues;
	file128.open("seq128.bin", std::ios::out | std::ios::binary);
	file_residues.open("residues.bin", std::ios::out | std::ios::binary);
	for (size_t i = 0; i < 3; i++)
	{
		file128.write((char*)&length128[i], sizeof(int));
	}
	offset = 0;
	unsigned int size128 = 0;
	unsigned int size128total = 0;
	unsigned int size128compr = 0;
	unsigned int size128comprtotal = 0;
	std::vector<std::pair<uint16_t, uint16_t>> rleVectors[3];
	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < compressVector[i].size(); j++)
		{
			std::pair<uint16_t, uint16_t> tempVal;
			if (compressVector[i].at(j).first == 128 && compressVector[i].at(j).second > 4)
			{
				size128 += compressVector[i].at(j).second;
				if (offset > 65535)
				{
					unsigned short temp = 65535;
					unsigned short temp_data = 0;
					while (offset > 65535) {
						file128.write((char*)&temp, sizeof(unsigned short));
						file128.write((char*)&temp_data, sizeof(unsigned short));
						tempVal.first = temp;
						tempVal.second = temp_data;
						rleVectors[i].push_back(tempVal);
						text_compress << std::dec << temp << " - " << temp_data << "\n";
						size128compr += 4;
						offset -= 65535;
					}
					temp = (unsigned short)offset;
					file128.write((char*)&temp, sizeof(unsigned short));
					file128.write((char*)&compressVector[i].at(j).second, sizeof(unsigned short));
					tempVal.first = temp;
					tempVal.second = compressVector[i].at(j).second;
					rleVectors[i].push_back(tempVal);
					size128compr += 4;
					text_compress << std::dec << temp << " - " << compressVector[i].at(j).second << "\n";
					offset = 0;
				}
				else 
				{
					unsigned short temp = (unsigned short)offset;
					file128.write((char*)&temp, sizeof(unsigned short));
					file128.write((char*)&compressVector[i].at(j).second, sizeof(unsigned short));
					tempVal.first = temp;
					tempVal.second = compressVector[i].at(j).second;
					rleVectors[i].push_back(tempVal);
					size128compr += 4;
					text_compress << std::dec << temp << " - " << compressVector[i].at(j).second << "\n";
					offset = 0;
				}
			}
			else
			{
				offset += compressVector[i].at(j).second;
			}
		}
		offset = 0;
		text_compress << "Size of 128 sequences is : " << size128 << " bytes; compressed to " << size128compr << "; ratio - " << (double)size128/size128compr << std::endl;
		size128total += size128;
		size128comprtotal += size128compr;
		size128 = 0;
		size128compr = 0;
	}

	for (size_t i = 0; i < 3; i++)
	{
		file128.write((char*)&lengthnon128[i], sizeof(int));
		//std::cout << lengthnon128[i] << std::endl;
		file_residues.write((char*)&lengthnon128[i], sizeof(int));
	}
	// store residues
	compressResFile(file128, frequencyTable, compressVector);

	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < compressVector[i].size(); j++)
		{
			if (compressVector[i].at(j).first != 128)
			{
				for (size_t k = 0; k < compressVector[i].at(j).second; k++)
				{
					file_residues.write((char*)&compressVector[i].at(j).first, sizeof(unsigned char));
				}
			}
			else if (compressVector[i].at(j).first == 128 && compressVector[i].at(j).second <= 4)
			{
				for (size_t k = 0; k < compressVector[i].at(j).second; k++)
				{
					file_residues.write((char*)&compressVector[i].at(j).first, sizeof(unsigned char));
				}
			}
		}
	}

	file_residues.flush();
	file_residues.close();

	file128.flush();
	file128.close();

	text_compress << "Total size of 128 sequences is : " << size128total << " bytes; compressed to " << size128comprtotal << "; ratio - " << (double)size128total / size128comprtotal << std::endl;
	
	text_compress.flush();
	text_compress.close();

	// compress all sequences and residues into a single file as a single stream
	compressAsStream(numLines, frequencyTable, rleVectors, compressVector);
	compressAsStreamElement("element0.bin", numLines, frequencyTable, (rleVectors + 0), (compressVector + 0));
	compressAsStreamElement("element1.bin", numLines, frequencyTable, (rleVectors + 1), (compressVector + 1));
	compressAsStreamElement("element2.bin", numLines, frequencyTable, (rleVectors + 2), (compressVector + 2));
	
	std::vector<std::pair<unsigned short, unsigned short>> compressRestore[3];
	
	std::vector<unsigned char> residuesVectors[3];
	
	// restore data structures to reproduce the original file
	decompressAsStream(numLines, compressRestore, residuesVectors);

	// restore file from sequences
	bool procSequence[3];
	int seqMarker[3] = { 0, 0, 0 };
	int resMarker[3] = { 0, 0, 0 };
	int curValue[3] = { 0, 0, 0 };
	int seqLength[3] = { 0, 0, 0 };
	for (size_t i = 0; i < 3; i++)
	{
		if (compressRestore[i].at(0).first == 0)
		{
			procSequence[i] = true;
			seqLength[i] = compressRestore[i].at(0).second;
		}
		else
		{
			procSequence[i] = false;
			seqLength[i] = compressRestore[i].at(0).first;
		}
	}
	std::ofstream frestored;
	frestored.open("fresy_restored.bin", std::ios::binary);
	for (size_t i = 0; i < numLines; i++)
	{
		for (size_t j = 0; j < 3; j++)
		{
			if (procSequence[j])
			{
				curValue[j] = 128;
				seqLength[j]--;
				if (seqLength[j] == 0)
				{
					seqMarker[j]++;
					if (seqMarker[j] == compressRestore[j].size()-1 && compressRestore[j].at(seqMarker[j]).second != 0)
					{
						procSequence[j] = true;
					}
					else if (seqMarker[j] == compressRestore[j].size() - 1 && compressRestore[j].at(seqMarker[j]).second == 0) {
						procSequence[j] = false;
					}
					if (seqMarker[j] == compressRestore[j].size()) continue;
					if (compressRestore[j].at(seqMarker[j]).first != 0)
					{
						seqLength[j] = compressRestore[j].at(seqMarker[j]).first;
						procSequence[j] = false;
					}
					else
					{
						seqLength[j] = compressRestore[j].at(seqMarker[j]).second;
						procSequence[j] = true;
					}
				}
			}
			else
			{
				if (resMarker[j] == residuesVectors[j].size()) continue;
				curValue[j] = residuesVectors[j].at(resMarker[j]);
				resMarker[j]++;
				seqLength[j]--;
				if (seqLength[j] == 0)
				{
					if (compressRestore[j].at(seqMarker[j]).second != 0)
					{
						procSequence[j] = true;
						seqLength[j] = compressRestore[j].at(seqMarker[j]).second;
					}
					else
					{
						seqMarker[j]++;
						procSequence[j] = false;
						if (seqMarker[j] >= compressRestore[j].size()) continue;
						seqLength[j] = compressRestore[j].at(seqMarker[j]).first;
					}
				}
			}
		}
		//frestored << std::hex << std::uppercase << curValue[0] << " " << curValue[1] << " " << curValue[2] << "\n";
		for (int k = 0; k < 3; k++) {
			frestored.write((char*)&curValue[k], sizeof(unsigned char));
		}
		//frestored.flush();
	}
	frestored.flush();
	frestored.close();
    std::cout << "Finished!\n";
	long before = get_file_size(argv[1]);
	long after = get_file_size("compressed.bin");
	std::cout << "File size before: " << before << "; size after: " << after << "; ratio: " << (double)before / after << std::endl;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
