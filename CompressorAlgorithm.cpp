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

void writeSymbolToCoder(HuffmanEncoder& enc, uint32_t symb) {
	int symbol = symb;
	if (symbol == EOF)
		return;
	if (symbol < 0 || symbol > 255)
		throw std::logic_error("Assertion error");
	enc.write(static_cast<uint32_t>(symbol));
}

void writeSequence(HuffmanEncoder& enc, std::pair<uint16_t, uint16_t>& pair) {
	uint8_t first = pair.first >> 8;
	uint8_t second = pair.first & 0xff;
	uint8_t third = pair.second >> 8;
	uint8_t fourth = pair.second & 0xff;
	uint8_t seqArr[4] = { first, second, third, fourth };
	for (size_t j = 0; j < 4; j++)
	{
		writeSymbolToCoder(enc, seqArr[j]);
	}
}

void writeSequence(HuffmanEncoder& enc, uint16_t f, uint16_t s) {
	uint8_t first = f >> 8;
	uint8_t second = f & 0xff;
	uint8_t third = s >> 8;
	uint8_t fourth = s & 0xff;
	uint8_t seqArr[4] = { first, second, third, fourth };
	for (size_t j = 0; j < 4; j++)
	{
		writeSymbolToCoder(enc, seqArr[j]);
	}
}

void compressAsStream(int numLines, std::vector<std::pair<uint16_t, uint16_t>>* rleVectors, std::vector<std::pair<int, int>>* compressVectors) {
	std::ofstream output;
	output.open("compressed.bin", std::ios::out | std::ios::binary);
	FrequencyTable freqs(std::vector<uint32_t>(256, 0));
	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < compressVectors[i].size(); j++)
		{
			if (compressVectors[i].at(j).first != 128 || (compressVectors[i].at(j).first == 128 && compressVectors[i].at(j).second <= 4))
				freqs.set(compressVectors[i].at(j).first, freqs.get(compressVectors[i].at(j).first) + compressVectors[i].at(j).second);
		}
	}
	for (size_t i = 0; i < 3; i++)
	{
		for (size_t j = 0; j < rleVectors[i].size(); j++)
		{
			uint8_t first = rleVectors[i].at(j).first >> 8;
			uint8_t second = rleVectors[i].at(j).first & 0xff;
			uint8_t third = rleVectors[i].at(j).second >> 8;
			uint8_t fourth = rleVectors[i].at(j).second & 0xff;
			freqs.set(first, freqs.get(first) + 1);
			freqs.set(second, freqs.get(second) + 1);
			freqs.set(third, freqs.get(third) + 1);
			freqs.set(fourth, freqs.get(fourth) + 1);
		}
	}
	CodeTree code = freqs.buildCodeTree();
	const CanonicalCode canonCode(code, freqs.getSymbolLimit());
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
			int curPos = 0;
			int resRun = 0;
			int next128 = 0;
			int posInVector = 0;
			// write first 128 sequence at the start
			writeSequence(enc, rleVectors[i].at(0));
			
			resRun = rleVectors[i].at(0).first;
			
			curPos += rleVectors[i].at(0).second;
			while (curPos != numLines)
			{
				// write residues
				while (resRun != 0)
				{
					if (compressVectors[i].at(posInVector).first != 128 || (compressVectors[i].at(posInVector).first == 128 
																			&& compressVectors[i].at(posInVector).second <= 4)) {
						for (size_t k = 0; k < compressVectors[i].at(posInVector).second; k++)
						{
							writeSymbolToCoder(enc, compressVectors[i].at(posInVector).first);
						}
						curPos += compressVectors[i].at(posInVector).second;
						resRun -= compressVectors[i].at(posInVector).second;
					}
					posInVector++;
					if (posInVector >= compressVectors[i].size())
					{
						break;
					}
				}
				// write next 128 sequence
				next128++;
				if (next128 == rleVectors[i].size()) break;
				if (next128 < rleVectors[i].size())
				{
					writeSequence(enc, rleVectors[i].at(next128));
					
					curPos += rleVectors[i].at(next128).second;
					if (next128 == rleVectors[i].size())
					{
						resRun = numLines - curPos;
					}
					else
					{
						resRun = rleVectors[i].at(next128).first;
					}
					
				}
				
			}
			if (curPos != numLines)
			{
				resRun = numLines - curPos;
				writeSequence(enc, resRun, 0);
				
				while (resRun != 0)
				{
					if (compressVectors[i].at(posInVector).first != 128 || (compressVectors[i].at(posInVector).first == 128
						&& compressVectors[i].at(posInVector).second <= 4)) {
						for (size_t k = 0; k < compressVectors[i].at(posInVector).second; k++)
						{
							writeSymbolToCoder(enc, compressVectors[i].at(posInVector).first);
						}
						curPos += compressVectors[i].at(posInVector).second;
						resRun -= compressVectors[i].at(posInVector).second;
					}
					posInVector++;
					if (posInVector >= compressVectors[i].size())
					{
						break;
					}
				}
			}
		}
		bout.finish();
	}
	catch (const char* msg) {
		std::cerr << msg << std::endl;
	}
}

uint8_t readSymbolFromDecoder(HuffmanDecoder& dec) {
	uint32_t symbol = dec.read();
	int b = static_cast<int>(symbol);
	if (std::numeric_limits<char>::is_signed)
		b -= (b >> 7) << 8;
	return static_cast<uint8_t>(b);
}

void decompressAsStream(int numLines, std::vector<std::pair<uint16_t, uint16_t>>* rleVectors, std::vector<unsigned char>* resVector) {
	std::ifstream input;
	input.open("compressed.bin", std::ios::in | std::ios::binary);
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
			int curPos = 0;
			int resRun = 0;
			int next128 = 0;
			std::pair<uint16_t, uint16_t> tempPair;
			uint8_t first = readSymbolFromDecoder(dec);
			uint8_t second = readSymbolFromDecoder(dec);
			uint8_t third = readSymbolFromDecoder(dec);
			uint8_t fourth = readSymbolFromDecoder(dec);
			tempPair.first = (first << 8) + second;
			tempPair.second = (third << 8) + fourth;
			rleVectors[i].push_back(tempPair);
			resRun = rleVectors[i].at(0).first;
			curPos += rleVectors[i].at(0).second;
			uint8_t temp;
			while (curPos != (numLines))
			{
				// process residues
				while (resRun != 0)
				{
					temp = readSymbolFromDecoder(dec);
					resVector[i].push_back(temp);
					resRun--;
					curPos++;
				}
				if (curPos == (numLines)) break;
				// process sequences
				uint8_t first = readSymbolFromDecoder(dec);
				uint8_t second = readSymbolFromDecoder(dec);
				uint8_t third = readSymbolFromDecoder(dec);
				uint8_t fourth = readSymbolFromDecoder(dec);
				tempPair.first = (first << 8) + second;
				tempPair.second = (third << 8) + fourth;
				rleVectors[i].push_back(tempPair);
				resRun = rleVectors[i].at(rleVectors[i].size()-1).first;
				curPos += rleVectors[i].at(rleVectors[i].size() - 1).second;
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
	for (size_t i = 0; i < 256; i++)
	{
		std::cout << i << "-" << frequencyTable.get(i) << std::endl;
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
		//file_residues << lengthnon128[i];
	}
	// store residues
	compressResFile(file128, frequencyTable, compressVector);

	file128.flush();
	file128.close();

	text_compress << "Total size of 128 sequences is : " << size128total << " bytes; compressed to " << size128comprtotal << "; ratio - " << (double)size128total / size128comprtotal << std::endl;
	
	text_compress.flush();
	text_compress.close();

	// compress all sequences and residues into a single file as a single stream
	compressAsStream(numLines, rleVectors, compressVector);
	
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
					if (seqMarker[j] == compressRestore[j].size()-1)
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
