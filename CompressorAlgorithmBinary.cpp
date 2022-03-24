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

void compressHuffman(char* outfilename) {
	std::ifstream in("residues.bin", std::ios::binary);
	FrequencyTable freqs(std::vector<uint32_t>(257, 0));
	while (true) {
		int b = in.get();
		if (b == EOF)
			break;
		if (b < 0 || b > 255)
			throw std::logic_error("Assertion error");
		freqs.increment(static_cast<uint32_t>(b));
	}
	freqs.increment(256);  // EOF symbol gets a frequency of 1
	CodeTree code = freqs.buildCodeTree();
	const CanonicalCode canonCode(code, freqs.getSymbolLimit());
	// Replace code tree with canonical one. For each symbol,
	// the code value may change but the code length stays the same.
	code = canonCode.toCodeTree();

	// Read input file again, compress with Huffman coding, and write output file
	in.clear();
	in.seekg(0);
	std::ofstream out(outfilename, std::ios::binary);
	BitOutputStream bout(out);
	try {

		// Write code length table
		for (uint32_t i = 0; i < canonCode.getSymbolLimit(); i++) {
			uint32_t val = canonCode.getCodeLength(i);
			// For this file format, we only support codes up to 255 bits long
			if (val >= 256)
				throw std::domain_error("The code for a symbol is too long");
			// Write value as 8 bits in big endian
			for (int j = 7; j >= 0; j--)
				bout.write((val >> j) & 1);
		}

		HuffmanEncoder enc(bout);
		enc.codeTree = &code;
		while (true) {
			// Read and encode one byte
			int symbol = in.get();
			if (symbol == EOF)
				break;
			if (symbol < 0 || symbol > 255)
				throw std::logic_error("Assertion error");
			enc.write(static_cast<uint32_t>(symbol));
		}
		enc.write(256);  // EOF
		bout.finish();
	}
	catch (const char* msg) {
		std::cerr << msg << std::endl;
	}
}

void decompressHuffman(char* outfilename) {
	std::ifstream in("residuesHuffman.bin", std::ios::binary);
	std::ofstream out(outfilename, std::ios::binary);
	BitInputStream bin(in);
	try {

		// Read code length table
		std::vector<uint32_t> codeLengths;
		for (int i = 0; i < 257; i++) {
			// For this file format, we read 8 bits in big endian
			uint32_t val = 0;
			for (int j = 0; j < 8; j++)
				val = (val << 1) | bin.readNoEof();
			codeLengths.push_back(val);
		}
		const CanonicalCode canonCode(codeLengths);
		const CodeTree code = canonCode.toCodeTree();

		HuffmanDecoder dec(bin);
		dec.codeTree = &code;
		while (true) {
			uint32_t symbol = dec.read();
			if (symbol == 256)  // EOF symbol
				break;
			int b = static_cast<int>(symbol);
			if (std::numeric_limits<char>::is_signed)
				b -= (b >> 7) << 8;
			out.put(static_cast<char>(b));
		}
	}
	catch (const char* msg) {
		std::cerr << msg << std::endl;
	}
}

void compressResFile(std::ofstream& output, FrequencyTable& ft, std::vector<std::pair<int, int>>* vector) {
	CodeTree code = ft.buildCodeTree();
	const CanonicalCode canonCode(code, ft.getSymbolLimit());
	// Replace code tree with canonical one. For each symbol,
	// the code value may change but the code length stays the same.
	code = canonCode.toCodeTree();

	BitOutputStream bout(output);
	try {

		// Write code length table
		for (uint32_t i = 0; i < canonCode.getSymbolLimit(); i++) {
			uint32_t val = canonCode.getCodeLength(i);
			// For this file format, we only support codes up to 255 bits long
			if (val >= 256)
				throw std::domain_error("The code for a symbol is too long");
			// Write value as 8 bits in big endian
			for (int j = 7; j >= 0; j--)
				bout.write((val >> j) & 1);
		}

		HuffmanEncoder enc(bout);
		enc.codeTree = &code;
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < vector[i].size(); j++)
			{
				if (vector[i].at(j).first != 128)
				{
					for (size_t k = 0; k < vector[i].at(j).second; k++)
					{
						enc.write(vector[i].at(j).first);
					}
				}
				else if (vector[i].at(j).first == 128 && vector[i].at(j).second <= 4)
				{
					for (size_t k = 0; k < vector[i].at(j).second; k++)
					{
						enc.write(vector[i].at(j).first);
					}
				}
			}
		}
		//enc.write(256);  // EOF
		bout.finish();
	}
	catch (const char* msg) {
		std::cerr << msg << std::endl;
	}
}

void decompressResFile(std::ifstream& input, std::vector<unsigned char> * resVector, int * len_res) {
	BitInputStream bin(input);
	try {

		// Read code length table
		std::vector<uint32_t> codeLengths;
		for (int i = 0; i < 256; i++) {
			// For this file format, we read 8 bits in big endian
			uint32_t val = 0;
			for (int j = 0; j < 8; j++)
				val = (val << 1) | bin.readNoEof();
			codeLengths.push_back(val);
		}
		const CanonicalCode canonCode(codeLengths);
		const CodeTree code = canonCode.toCodeTree();

		HuffmanDecoder dec(bin);
		dec.codeTree = &code;
		for (size_t i = 0; i < 3; i++)
		{
			for (size_t j = 0; j < len_res[i]; j++)
			{
				uint32_t symbol = dec.read();
				unsigned char temp;
				int b = static_cast<int>(symbol);
				if (std::numeric_limits<char>::is_signed)
					b -= (b >> 7) << 8;
				resVector[i].push_back(static_cast<unsigned char>(b));
			}
		}
	}
	catch (const char* msg) {
		std::cerr << msg << std::endl;
	}
}

long get_file_size(const char* filename) {
	struct stat statv;
	int res = stat(filename, &statv);
	return res == 0 ? statv.st_size : -1;
}

int main(int argc, char * argv[])
{
    std::ifstream infile;
	if (argc != 2)
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
	// calculate total length of non-compressed sequences in each column
	FrequencyTable frequencyTable(std::vector<uint32_t>(256, 0));
	unsigned int offset = 0;
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
	// store 128-based sequences
	std::ofstream text_compress("compress.txt", std::ios::out);
	std::ofstream file128;
	file128.open("seq128.bin", std::ios::out | std::ios::binary);
	for (size_t i = 0; i < 3; i++)
	{
		file128.write((char*)&length128[i], sizeof(int));
	}
	offset = 0;
	unsigned int size128 = 0;
	unsigned int size128total = 0;
	unsigned int size128compr = 0;
	unsigned int size128comprtotal = 0;
	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < compressVector[i].size(); j++)
		{
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
						text_compress << std::dec << temp << " - " << temp_data << "\n";
						size128compr += 4;
						offset -= 65535;
					}
					temp = (unsigned short)offset;
					file128.write((char*)&temp, sizeof(unsigned short));
					file128.write((char*)&compressVector[i].at(j).second, sizeof(unsigned short));
					size128compr += 4;
					text_compress << std::dec << temp << " - " << compressVector[i].at(j).second << "\n";
					offset = 0;
				}
				else 
				{
					unsigned short temp = (unsigned short)offset;
					file128.write((char*)&temp, sizeof(unsigned short));
					file128.write((char*)&compressVector[i].at(j).second, sizeof(unsigned short));
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
		//file_residues << lengthnon128[i];
	}
	// store residues
	compressResFile(file128, frequencyTable, compressVector);

	file128.flush();
	file128.close();

	text_compress << "Total size of 128 sequences is : " << size128total << " bytes; compressed to " << size128comprtotal << "; ratio - " << (double)size128total / size128comprtotal << std::endl;
	
	text_compress.flush();
	text_compress.close();
	// restore all sequences
	std::ifstream in128;
	in128.open("seq128.bin", std::ios::in | std::ios::binary);
	int len128[3] = {0, 0, 0};
	for (size_t i = 0; i < 3; i++)
	{
		in128.read((char*)&len128[i], sizeof(int));
	}
	std::vector<std::pair<unsigned short, unsigned short>> compressRestore[3];
	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < len128[i]; j++)
		{
			std::pair<unsigned short, unsigned short> tempPair;
			unsigned short temp;
			in128.read((char*)&temp, sizeof(unsigned short));
			tempPair.first = temp;
			in128.read((char*)&temp, sizeof(unsigned short));
			tempPair.second = temp;
			compressRestore[i].push_back(tempPair);
		}
	}
	
	// restore residues
	int len_res[3] = { 0, 0, 0 };
	for (size_t i = 0; i < 3; i++)
	{
		in128.read((char*)&len_res[i], sizeof(int));
	//	//std::cout << len_res[i] << std::endl;
	}
	std::vector<unsigned char> residuesVectors[3];
	
	decompressResFile(in128, residuesVectors, len_res);
	in128.close();
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
					if (seqMarker[j] == compressRestore[j].size())
					{
						procSequence[j] = false;
						continue;
					}
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
						seqLength[j] = compressRestore[j].at(seqMarker[j]).first;
					}
				}
			}
		}
		//frestored << std::hex << std::uppercase << curValue[0] << " " << curValue[1] << " " << curValue[2] << "\n";
		for (int k = 0; k < 3; k++)
			frestored.write((char*)&curValue[k], sizeof(unsigned char));
	}
	frestored.flush();
	frestored.close();
    std::cout << "Finished!\n";
	long before = get_file_size(argv[1]);
	long after = get_file_size("seq128.bin");
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
