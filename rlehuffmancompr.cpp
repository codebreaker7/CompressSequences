#include "rlehuffmancompr.h"
#include "CanonicalCode.hpp"

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

void decompressResFile(std::ifstream& input, std::vector<unsigned char>* resVector, int* len_res) {
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

void compressAsStream(int numLines, FrequencyTable& freqs, std::vector<std::pair<uint16_t, uint16_t>>* rleVectors, std::vector<std::pair<int, int>>* compressVectors) {
	std::ofstream output;
	output.open("compressed.bin", std::ios::out | std::ios::binary);
	
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
	std::vector<char> sequence;
	// write corresponding symbols as a code
	for (size_t i = 0; i < 256; i++)
	{
		sequence = code.getCode(i);
		std::cout << i << "-";
		for (size_t j = 0; j < sequence.size(); j++)
		{
			if (sequence.at(j) == 0)
			{
				std::cout << "0";
			}
			else
			{
				std::cout << "1";
			}
		}
		std::cout << std::endl;
	}
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

void compressAsStreamElement(const char* outfile, int numLines, FrequencyTable& freqs, std::vector<std::pair<uint16_t, uint16_t>>* rleVectors, std::vector<std::pair<int, int>>* compressVectors) {
	std::ofstream output;
	output.open(outfile, std::ios::out | std::ios::binary);

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

		for (size_t i = 0; i < 1; i++) // treat this parameter as an array
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
				resRun = rleVectors[i].at(rleVectors[i].size() - 1).first;
				curPos += rleVectors[i].at(rleVectors[i].size() - 1).second;
			}
		}
	}
	catch (const char* msg) {
		std::cerr << msg << std::endl;
	}
}

FrequencyTable loadFrequencyTable(const char* filename) {
	FrequencyTable ft(std::vector<uint32_t>(256, 0));
	std::ifstream freq_file;
	freq_file.open(filename, std::ios::in);
	uint32_t val, freq;
	while (freq_file >> val >> freq)
	{
		if (freq < 5) freq = 5;
		ft.set(val, freq);
	}
	freq_file.close();
	return ft;
}