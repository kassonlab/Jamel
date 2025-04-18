#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
// /home/vs_cpp venv  source /home/vs_cpp/bin/activate
using namespace std;
//  g++ ./Chimeragenesis/Scripts/test.cpp -o test
std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<std::map<std::string,std::map<std::tuple<long long,long long>, std::tuple<long long,long long>>>>> single_cut_chimera_generator(std::string base_label,std::string base_seq, std::string partner_label,std::string partner_seq,float lower_bound,float upper_bound) {
    std::size_t base_len=base_seq.length();
    std::vector<int> base_cuts;
    std::vector<std::string> chi_labels;
    std::vector<std::string> chi_seqs;
    std::vector<std::map<std::string,std::map<std::tuple<long long,long long>, std::tuple<long long,long long>>>> inheritances;

    for (int i = std::round(lower_bound*base_len); i < std::round(upper_bound*base_len); ++i) {
        base_cuts.push_back(i);
    }
    std::size_t partner_len=partner_seq.length();
    std::vector<int> partner_cuts;

    for (int i = std::round(lower_bound*partner_len); i < std::round(upper_bound*partner_len); ++i) {
        partner_cuts.push_back(i);
    }
    for (auto base_cut : base_cuts) {
        for (auto partner_cut: partner_cuts) {
            std::string base_splice=base_seq.substr(0,base_cut);
            std::string partner_splice=partner_seq.substr(partner_cut);
            std::string chi_seq=base_splice+partner_splice;
            long long splice_len=base_splice.length();
            long long chimera_len=chi_seq.length();

            std::map<std::string,std::map<std::tuple<long long,long long>, std::tuple<long long,long long>>> inheritance;
            std::tuple<long long,long long> base_key=std::make_tuple(0,splice_len);
            std::tuple<long long,long long> base_value=std::make_tuple(0,splice_len);
            std::map<std::tuple<long long,long long>, std::tuple<long long,long long>> base_map;
            base_map[base_key]=base_value;
            inheritance[base_label]=base_map;

            std::tuple<long long,long long> partner_key=std::make_tuple(partner_cut,partner_len);
            std::tuple<long long,long long> partner_value=std::make_tuple(splice_len,chimera_len);
            std::map<std::tuple<long long,long long>, std::tuple<long long,long long>> partner_map;
            partner_map[partner_key]=partner_value;
            inheritance[partner_label]=partner_map;
            
            chi_labels.push_back(base_label+"_0_"+std::to_string(base_cut)+"_"+partner_label+"_"+std::to_string(partner_cut)+"_"+std::to_string(partner_len));
            chi_seqs.push_back(chi_seq);
            inheritances.push_back(inheritance);
    }
}
    return std::make_tuple(chi_labels,chi_seqs,inheritances);
}

PYBIND11_MODULE(chi_cpp, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("single_cut_chimera_generator", &single_cut_chimera_generator, "Creates all possible combinations between a base protein at the n-terminus and a partner at the c-term");
}