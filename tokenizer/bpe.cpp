
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <tuple>
#include <utility>
#include <numeric>
#include <cmath>
#include <map>
#include <unordered_map>
#include <fstream>
#include "nlohmann/json.hpp"
#include "tokens_model.hpp"
#include <mutex>
#include <thread>
#include <future>
#include <execution>
#include <atomic>
#include <algorithm>
#include <iterator>
#include <utility>


using json = nlohmann::json;
typedef std::uint16_t TokenType;
typedef std::tuple<TokenType, TokenType> kmer;

const uint N_HELP_TOKENS = 6;
const uint MAX_N_TOKENS = 65535;



std::map<std::string, TokenType> alphabet = {
    {"[UNK]", 0},
    {"[CLS]", 1},
    {"[SEP]", 2},
    {"[PAD]", 3},
    {"[MASK]", 4},
    {"~", 5},
    {"N", 6},
    {"A", 7},
    {"C", 8},
    {"G", 9},
    {"T", 10}
};


template<typename T>
struct tuple_compare {
    bool operator()(const std::tuple<T, T>& a, const std::tuple<T, T>& b) const {
        if (std::get<0>(a) == std::get<0>(b)) {
            return std::get<1>(a) < std::get<1>(b);
        }
        return std::get<0>(a) < std::get<0>(b);
    }
};

// Custom hash function for std::tuple
template<typename T>
struct tuple_hash {
    template<typename Tuple, TokenType N>
    struct hasher {
        static void hash(Tuple const& t, TokenType& seed) {
            hasher<Tuple, N - 1>::hash(t, seed);
            seed ^= std::hash<typename std::tuple_element<N - 1, Tuple>::type>()(std::get<N - 1>(t)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
    };

    template<typename Tuple>
    struct hasher<Tuple, 1> {
        static void hash(Tuple const& t, TokenType& seed) {
            seed ^= std::hash<typename std::tuple_element<0, Tuple>::type>()(std::get<0>(t)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
    };

    template<typename... Args>
    TokenType operator()(std::tuple<Args...> const& t) const {
        TokenType seed = 0;
        hasher<std::tuple<Args...>, sizeof...(Args)>::hash(t, seed);
        return seed;
    }
};


void get_sequences_trf(const std::string& trf_file_name, std::vector<std::string>& seqs) {
    std::ifstream fh(trf_file_name);
    if (fh.is_open()) {
        std::string line;
        while (std::getline(fh, line)) {
            std::vector<std::string> data;
            size_t pos = 0;
            while ((pos = line.find('\t')) != std::string::npos) {
                data.push_back(line.substr(0, pos));
                line.erase(0, pos + 1);
            }
            seqs.push_back(data[14]);
        }
        fh.close();
    }
}

void get_sequences_reads(const std::string& reads_file_name, std::vector<std::string>& seqs) {
    std::ifstream fh(reads_file_name);
    if (fh.is_open()) {
        std::string line;
        while (std::getline(fh, line)) {
            seqs.push_back(line);
        }
        fh.close();
    }
}

std::string get_dataset(const std::vector<std::string>& seqs) {
    std::stringstream ss;
    bool first = true;
    for (const auto& s : seqs) {
        if (!first) {
            ss << "~";
        }
        first = false;
        ss << s;
    }
    return ss.str();
}

// in our case maximal tokens for std::uint16_t is 65535
std::vector<TokenType> convert_to_vector(std::string& dataset) {
    std::vector<TokenType> seq;
    seq.reserve(dataset.size()); // Reserve space
    for (auto x : dataset) {
        if (x == '\n') {
            x = '~';
        }
        if (alphabet.find(std::string(1, x)) != alphabet.end()) {
            seq.push_back(alphabet.at(std::string(1, x)));
        } else {
            seq.push_back(alphabet.at("[UNK]"));
        }
    }
    return seq;
}


std::string token_type_to_string(TokenType token, const std::map<std::string, TokenType>& alphabet, const std::map<TokenType, kmer>& tokens) {
    
    // Check if the token is in the initial alphabet
    for (const auto& element : alphabet) {
        if (element.second == token) {
            return element.first;
        }
    }

    // Check if the token is in the tokens map
    auto token_iter = tokens.find(token);
    if (token_iter != tokens.end()) {
        kmer kmer_ = token_iter->second;
        TokenType token1 = std::get<0>(kmer_);
        TokenType token2 = std::get<1>(kmer_);
        return token_type_to_string(token1, alphabet, tokens) + token_type_to_string(token2, alphabet, tokens);
    }

    // If the token is not found in either map, return an empty string
    return "";
}


void compute_freqs(const std::vector<TokenType>& seq, std::vector<std::atomic_size_t>& c, std::vector<std::thread>& threads, size_t n_threads, TokenType L) {

    std::cout << "Computing frequencies" << std::endl;

    auto worker = [&](size_t start, size_t end) {
        for (size_t i = start; i < end; ++i) {
            if (seq[i] > N_HELP_TOKENS && seq[i + 1] > N_HELP_TOKENS) {
                ++c[seq[i] * L + seq[i + 1]];
            }
        }
    };
    size_t chunk_size = (seq.size() - 1) / n_threads;

    for (size_t i = 0; i < n_threads; ++i) {
        size_t start = i * chunk_size;
        size_t end = (i == n_threads - 1) ? seq.size() - 1 : start + chunk_size;
        threads.emplace_back(worker, start, end);
    }

    for (auto &t : threads) {
        t.join();
    }
}

std::pair<std::size_t, kmer> found_max(const std::vector<std::atomic_size_t>& c, TokenType L) {
    std::cout << "find max " << std::endl;
    size_t max_count = 0;
    kmer rep;
    for (size_t i=0; i < c.size(); ++i) {
        size_t current_count = c[i].load();
        if (current_count > max_count) {
            max_count = current_count;
            rep = std::make_tuple(i / L, i % L);
        }
    }
    return std::make_pair(max_count, rep);
}

void transform_data(std::vector<TokenType> &seq, std::vector<kmer> &merged, std::map<TokenType, kmer> &tokens, std::unordered_map<kmer, TokenType, tuple_hash<TokenType>> &rev_tokens, size_t max_tokens, std::vector<bool>& to_replace, kmer& rep, size_t tf, std::vector<TokenType>& new_seq,  TokenType L, uint k=2) {
    
    

    std::fill(to_replace.begin(), to_replace.end(), false);
    for (size_t i = 0; i < seq.size() - k + 1; i++) {
        if (seq[i] > N_HELP_TOKENS && seq[i+1] > N_HELP_TOKENS) {
            kmer kmer_ = std::make_tuple(seq[i], seq[i+1]);
            if (kmer_ == rep) {
                to_replace[i] = true;
            }
        }
    }
    
    new_seq.clear();
    new_seq.reserve(seq.size());
    for (size_t i = 0; i < seq.size();) {
        if (to_replace[i]) {
            new_seq.push_back(L);
            i += 2;
        } else {
            new_seq.push_back(seq[i]);
            i += 1;
        }
    }

    seq = std::vector<TokenType>(new_seq.begin(), new_seq.end());
    
}

int main(int argc, char* argv[]) {

    
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file> <output_model_file> <max_tokens> <threads>" << std::endl;
        return 1;
    }

    std::string file_name = argv[1];
    std::string output_file = argv[2];
    std::string output_model_file = argv[3];
    size_t max_tokens = std::stoul(argv[4]);
    size_t n_threads = std::stoul(argv[5]);

    if (max_tokens > MAX_N_TOKENS) {
        std::cout << "Max tokens must be less than 65535" << std::endl;
        return 1;
    }

    std::vector<std::string> seqs;
    std::cout << "read file" << std::endl;
    get_sequences_reads(file_name, seqs);
    std::cout << "get dataset" << std::endl;
    std::string dataset = get_dataset(seqs);
    seqs.clear();
    
    std::vector<TokenType> seq = convert_to_vector(dataset);
    dataset.clear();

    std::vector<kmer> merged;
    uint k = 2;
    TokenType L = alphabet.size();
    std::map<TokenType, kmer> tokens;
    std::unordered_map<kmer, TokenType, tuple_hash<TokenType>> rev_tokens;

    std::vector<TokenType> new_seq;
    std::vector<bool> to_replace(seq.size(), false);
    while (true) {
        std::cout << "Tokens " << L << " count reps ";

       

        // Allocate memory for the matrix
        std::vector<std::atomic_size_t> c(L * L);
        // Initialize the elements to zero
        for (auto &elem : c) {
            elem.store(0);
        }

        std::vector<std::thread> threads;
        compute_freqs(seq, c, threads, n_threads, L);

        auto max_result = found_max(c, L);
        size_t tf = max_result.first;
        kmer rep = max_result.second;
        if (tf == 1) {
            break;
        }

        merged.push_back(rep);
        tokens[L] = rep;
        rev_tokens[rep] = L;

        transform_data(seq, merged, tokens, rev_tokens, max_tokens, to_replace, rep, tf, new_seq, L, k);


        std::string token1 = token_type_to_string(std::get<0>(rep), alphabet, tokens);
        std::string token2 = token_type_to_string(std::get<1>(rep), alphabet, tokens);
        std::cout << token1 << " " << token2 << " " << tf << " : "<< "replace: " << seq.size() << " -> " << new_seq.size() << std::endl;

        L += 1;
        if (max_tokens && L > max_tokens) {
            break;
        }
        
    }
    // compute tokens as strings
    std::map<TokenType, std::string> tokens_str_map;
    for (const auto& element : tokens) {
        std::string token_string = token_type_to_string(element.first, alphabet, tokens);
        tokens_str_map[element.first] = token_string;
    }
    for (const auto& element : alphabet) {
        tokens_str_map[element.second] = element.first;
    }

    // std::ofstream out_file(output_file);
    // if (out_file.is_open()) {
    //     for (const auto& element : seq) {
    //         if (element == 5) {
    //             out_file << "\n";
    //         } else {
    //             out_file << tokens_str_map.at(element) << " ";
    //         }
    //     }
    //     out_file << std::endl;
    //     out_file.close();
    // }

    nlohmann::ordered_json json_data = get_json(tokens_str_map, tokens);

    std::ofstream configFile(output_model_file);
    configFile << std::setw(2) << json_data << std::endl;
    configFile.close();

    std::cout << "Config saved to config.json" << std::endl;

    return 0;
}
