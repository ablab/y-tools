#include "verify.hpp"
#include "permutation.hpp"

void Permutation::ReadFromFile(std::string filename) {
    std::ifstream perm_fhandler(filename.c_str());
    VERIFY_MSG(!perm_fhandler.fail(), "Permutation file " << filename << " was not found");
    size_t index1 = 0;
    size_t index2 = 0;
    while(!perm_fhandler.eof()) {
        std::string line;
        getline(perm_fhandler, line);
        if(line == "")
            break;

        std::stringstream ss(line);
        ss >> index2;

        direct_[index1] = index2;
        reverse_[index2] = index1;
        index1++;
    }
    perm_fhandler.close();
}

std::ostream& operator<<(std::ostream &out, const Permutation &permutation) {
    out << "Direct: ";
    for(size_t i = 0; i < permutation.Size(); i++)
        out << i << ": " << permutation.Direct()[i] << ", ";
    out << std::endl;
    out << "Reverse: ";
    for(size_t i = 0; i < permutation.Size(); i++)
        out << i << ": " << permutation.Reverse()[i] << ", ";
    return out;
}