#pragma once

#include <iostream>
#include <set>

inline char get_complementary(char nucl){
	if(nucl == 'A' || nucl == 'a')
		return 'T';
	if(nucl == 'T' || nucl == 't')
		return 'A';
	if(nucl == 'C' || nucl == 'c')
		return 'G';
	if(nucl == 'G' || nucl == 'g')
		return 'C';
	if(nucl == 'N' || nucl == 'n')
		return 'A';
    std::cout << "Char " << nucl << " is not from nucleotide alphabet" << std::endl;
	VERIFY(false);
	return 'A';
}

inline std::string reverse_complementary(std::string seq) {
	std::string rc_seq = seq;
	for(size_t i = 0; i < seq.size(); i++)
		rc_seq[seq.size() - i - 1] = get_complementary(seq[i]);
	return rc_seq;
}

inline size_t HammingDistance(std::string s1, std::string s2) {
	VERIFY(s1.size() == s2.size());
	size_t dist = 0;
	for(size_t i = 0; i < s1.size(); i++)
		if(s1[i] != s2[i])
			dist++;
	return dist;
}
