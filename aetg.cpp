#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <random>
#include <time.h>
#include <numeric>
#include <fstream>

using namespace std;

/*
    Function to generate random number within a range.
*/
int rand_num(int from, int to) {
    random_device dev;
    mt19937 rng(dev());
    uniform_int_distribution<mt19937::result_type> dist6(from, to);
    return dist6(rng);
}

/*
    Function to swap values between variables.
*/
void swap(int &x, int &y) {
    int temp = x;
    x = y;
    y = temp;
}

/*
    Hash implemenation for a pair. Allows use of a pair as a key in a hash table.
*/
struct Pair_Hash {
	size_t operator () (pair<int, int> const &pair) const {
		size_t first = hash<int>()(pair.first);
		size_t second = hash<int>()(pair.second);

		return (first ^ second) + second * second;
	}
};


class Pairs {
    private:
       unordered_set<pair<int, int>, Pair_Hash> pairs;          // essentially a hash table to store pairs
       unordered_map<int, unordered_set<int>> all_pairs;        // another hash table to store individual values and all potential values they can be paired with
       int uncovered_count;                                     // count of all uncovered pairs
       int total_pairs;                                         // count of all pairs before any have been covered

    public:
        Pairs(int factors, int levels, const vector<vector<int>>& table) {  //constructor that generates all possible pairs
            for(int i = 0; i < factors; i++) {                                  
                for(int j = 0; j < levels; j++) {
                    for(int k = i + 1; k < factors; k++) {
                        for(int l = 0; l < levels; l++) {
                            pairs.insert(pair<int, int>(i * levels + j, k * levels + l));
                        }
                    }       
                }
            }
            uncovered_count = pairs.size();
            total_pairs = pairs.size();

            int columns = table.size();
            for (int i = 0; i < columns; i++) {
                for (int j = 0; j < table.at(i).size(); j++) {
                    unordered_set<int> right_hand_side;
                    for (int k = 0; k < columns; k++) {
                        if (k != i) {
                            for (int l = 0; l < table.at(k).size(); l++) {
                                right_hand_side.insert(table.at(k).at(l));
                            }
                        }       
                    }
                    all_pairs.insert(pair<int, unordered_set<int>>(table.at(i).at(j), right_hand_side));
                }
            }
        }

        pair<int, int> getRandomPair() {            //vgenerate a random pair from the hash table
            auto it = begin(pairs);
            advance(it, rand_num(0, pairs.size() - 1));
            return *it;
        }

        bool isUncovered(pair<int, int> p) {                // checks if pairs is uncovered
            bool is_uncovered = pairs.find(p) != pairs.end();
            return is_uncovered;
        }

        void removePair(pair<int, int> p) {             // removes pair from both hash tables, decrements uncovered count
            if (pairs.erase(p) == 0) {
                return;
            }
            all_pairs[p.first].erase(p.second);
            all_pairs[p.second].erase(p.first);
            if (all_pairs[p.first].size() == 0) { 
                all_pairs.erase(p.first);
            }
            if (all_pairs[p.second].size() == 0) {
                all_pairs.erase(p.second);
            }
            uncovered_count--;
        }
       
        int getUncoveredPairCount() {               // get uncovered count
            return uncovered_count;
        }

        int getTotalPairs() {                       // get max pair count
            return total_pairs;
        }

        int findFirstBestFactor(const vector<int>& factors) {       // used when finding the first factor for a candidate
            int best_factor = -1, rand_choice, size;
            int best_factor_size = 0;
            for (auto x : factors) {
                size = all_pairs[x].size();   // assign pair count
                if (size > best_factor_size) {  // check if pair count is highest
                    best_factor_size = size;
                    best_factor = x;
                }
                else if (size == best_factor_size) {    // if tied, generate random number to break tie
                    rand_choice = rand_num(0, 1);
                    if (rand_choice == 0) {
                        best_factor = x;
                        best_factor_size = size;
                    }
                }
            }
            if (best_factor_size == 0) {    // in case no factor was found with >0 pairs
                    best_factor = factors.at(rand_num(0, factors.size() - 1));
            }
            return best_factor;
        }

        int findBestFactor(int index ,const vector<int>& candidate, const vector<int>& factors, const vector<vector<int>>& table, unordered_map<int, int>& map) {  // finds best factor, after at least the first factor has been determined
            int best_factor = -1, rand_choice, size;
            int uncovered, best_uncovered = 0;
            int left, right;
            
            for (auto x : factors) {    
                uncovered = 0;  // reset counter to zero for every factor
                for (auto y : candidate) {
                    if (x > y) {            // check to see which value of the pair is smaller
                        right = x;
                        left = y;
                    }
                    else {
                        right = y;
                        left = x;
                    }
        
                    if (isUncovered(pair<int, int>(left, right))) { // if pair is uncoverec, increment counter
                        uncovered++;
                    }
                }

                if (uncovered > best_uncovered) { // make sure the best factor is known
                    best_uncovered = uncovered;
                    best_factor = x;
                }
                else if (uncovered == best_uncovered) { // random number is generated to break ties
                    rand_choice = rand_num(0, 1);
                    if (rand_choice == 0) {
                        best_factor = x;
                        best_uncovered = uncovered;
                    }
                }
            }

                if (best_uncovered == 0) {
                    best_factor = factors.at(rand_num(0, factors.size() - 1));
                }
                return best_factor;
                
            }


            int checkUncovered(const vector<int>& candidate) { // check to see if a pair is uncovered
                int total = 0, left, right;
                for (int i = 0; i < candidate.size() - 1; i++) {
                    for (int j = i + 1; j < candidate.size(); j++) {
                        left = candidate.at(i);
                        right = candidate.at(j);
                        if (left > right) {
                            swap(left, right);
                        }
                        pair<int, int> p(left, right);
                        if (isUncovered(p)) {
                            total++;
                        }
                    }
                }
                 return total;
            }
           

        void removeCandidatePairs(const vector<int>& candidate) { // removes uncovered pairs from hash table
            int left, right;
            for (int i = 0; i < candidate.size() - 1; i++) {
                for (int j = i + 1; j < candidate.size(); j++) {
                    left = candidate.at(i);
                    right = candidate.at(j);
                    if (left > right) {
                        swap(left, right);
                    }
                    pair<int, int> p(left, right);
                    removePair(p);
                }
            }
        }        
};

unordered_map<int, int> generateComponentMap(int factors, int levels) { // generate a map, used to determing which set of factors a value belongs to
    unordered_map<int, int> map;
    int f = 0;
    for (int i = 0; i < factors; i++) {
        for (int j = 0; j < levels; j++) {
            map.insert(pair<int, int>(f, i));
            f++;
        }
    }
    return map;
}

vector<vector<int>> generateComponentTable(int factors, int levels) { // generates a table with all factors and levels
    vector<vector<int>> table;  
    int f = 0;
    for (int i = 0; i < factors; i++) {
        vector<int> column;
        for (int j = 0; j < levels; j++) {
            column.push_back(f);
            f++;
        }
        table.push_back(column);
    }
    return table;
}

int main()
{ 
    int factors, levels;

    cout << "Enter a value for the amount of factors: "; cin >> factors;
    cout << "Enter a value for the amount of levels: "; cin >> levels;

    unordered_map<int, int> component_map = generateComponentMap(factors, levels);
    vector<vector<int>> component_table = generateComponentTable(factors, levels);
    

    int suite_count = 0, candidate_count = 0, rand_choice;
    int total_pairs, average = 0, best, worst = 0;
    int candidate_cover_amount, best_candidate_cover_amount = 0;
    double total_time = 0;
    vector<int> candidate, best_candidate;
    vector<vector<int>> suite, best_suite;
    vector<int> factor_range(factors);                  // allocate size for vector
    iota(factor_range.begin(), factor_range.end(), 0); // fill vector with a range, starting at 0
    
    
    while (suite_count < 100) {             // generates 100 test suites
        clock_t begin = clock();            // begin timer
        suite.clear();
        Pairs pairs(factors, levels, component_table);  // create object to keep track of pairs
        total_pairs = pairs.getTotalPairs();
        while (pairs.getUncoveredPairCount() > 0) {     // loop until all pairs have been covered
           best_candidate.clear();
           best_candidate_cover_amount = 0;
           
            while (candidate_count < 50) {  // generate 50 candidates
                candidate.clear();  // reset candidate
                random_shuffle(factor_range.begin(), factor_range.end());   // randomly shuffle factors
                for (int i = 0; i < factor_range.size(); i++) { // determine factors
                        int best_factor;
                        if (i == 0) {
                            best_factor = pairs.findFirstBestFactor(component_table.at(factor_range.at(i))); // get best first factor
                            candidate.push_back(best_factor);
                        }
                        else {
                            best_factor = pairs.findBestFactor(factor_range.at(i), candidate, component_table.at(factor_range.at(i)), component_table, component_map); // get best factor after the first has been found
                            candidate.push_back(best_factor);
                        }
                }
            
                candidate_cover_amount = pairs.checkUncovered(candidate); // get the amount of pairs a candidate covers
                if (candidate_cover_amount > best_candidate_cover_amount) { // check to see if cover amount is best
                    best_candidate_cover_amount = candidate_cover_amount;
                    best_candidate = candidate;
                }
                else if (candidate_cover_amount == best_candidate_cover_amount) { // random number generator to break ties
                    rand_choice = rand_num(0, 1);
                    if (rand_choice == 0) {
                        best_candidate_cover_amount = candidate_cover_amount;
                        best_candidate = candidate;
                    }
                }
                candidate_count++;
                
            } 
            pairs.removeCandidatePairs(best_candidate);  // remove all pairs covered
            suite.push_back(best_candidate); // add candidate to suite
            candidate_count = 0; 
        }
        clock_t end = clock(); // end timer
        total_time += double(end - begin) / CLOCKS_PER_SEC; // calculate time spent finding suite
        if (suite_count == 0) { // set first initial suite
            best = suite.size();
            best_suite = suite;
        }
        else {
            if (suite.size() < best) { // check if suite is best
                best_suite = suite;
                best = best_suite.size();
            }
        }
        average += suite.size();
        if (suite.size() > worst){
            worst = suite.size();
        }
        
        suite_count++;
    }

    cout << "Best result: " << best << endl;
    cout << "Average result: " << average / 100 << endl;
    cout << "Worst result: " << worst << endl;
    cout << "Average Execution time: " << double(total_time / 100.0) << endl;

    ofstream fout;
    fout.open("best.txt");
    fout << best_suite.size() << endl << endl;
    for (auto x : best_suite) {
        sort(x.begin(), x.end());
        for (auto y : x) {
            fout << y << " ";
        }
        fout << endl;
    }
    fout.close();
    return 0; 
}