#pragma once

#include <iostream>
#include <fstream>
#include "FrequencyTable.hpp"
#include "HuffmanCoder.hpp"


void compressHuffman(char* outfilename);
void decompressHuffman(char* outfilename);
void compressResFile(std::ofstream& output, FrequencyTable& ft, std::vector<std::pair<int, int>>* vector);
void decompressResFile(std::ifstream& input, std::vector<unsigned char>* resVector, int* len_res);
void writeSymbolToCoder(HuffmanEncoder& enc, uint32_t symb);
void writeSequence(HuffmanEncoder& enc, std::pair<uint16_t, uint16_t>& pair);
void writeSequence(HuffmanEncoder& enc, uint16_t f, uint16_t s);
void compressAsStream(int numLines, std::vector<std::pair<uint16_t, uint16_t>>* rleVectors, std::vector<std::pair<int, int>>* compressVectors);
uint8_t readSymbolFromDecoder(HuffmanDecoder& dec);
void decompressAsStream(int numLines, std::vector<std::pair<uint16_t, uint16_t>>* rleVectors, std::vector<unsigned char>* resVector);
