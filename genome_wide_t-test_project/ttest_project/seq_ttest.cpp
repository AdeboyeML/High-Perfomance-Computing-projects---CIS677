#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator
#include <random>
#include <tuple>
#include <utility>
#include <cmath>
#include <numeric>      // std::iota
#include <algorithm> 
#include <chrono>

using namespace std;

template<typename T>
std::vector<T> flatten(const std::vector<std::vector<T>>& orig)
{
    std::vector<T> ret;
    for (const auto& v : orig)
        ret.insert(ret.end(), v.begin(), v.end());
    return ret;
}

template <typename T>
vector<size_t> sort_indexes(const vector<T>& v) {

    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using std::stable_sort instead of std::sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    stable_sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

    return idx;
}

vector<vector<string>> parseCsvToVector(string fileName) {
    std::ifstream csv(fileName);
    std::string line;
    std::vector <std::vector<std::string>> matRows;

    if (csv.is_open()) {
        for (std::string row_line; std::getline(csv, row_line);)
        {
            matRows.emplace_back();
            std::istringstream row_stream(row_line);
            for (std::string column; std::getline(row_stream, column, ',');)
                matRows.back().push_back(column);

        }
        csv.close();
    }
    else {
        cout << "Unable to open file";
    }

    /*for (size_t i = 0; i < matRows.size(); i++) {
        for (size_t j = 0; j < matRows[i].size(); j++) {
            cout << matRows[i][j] << " ";
        }
        cout << endl;
    }*/
    int Nrows = matRows.size();
    int Ncols = matRows[0].size();

    cout << "number of rows: " << Nrows << endl;
    cout << "number of columns: " << Ncols << endl;
    return matRows;

}

void write_csv(vector<vector<string>> dataset, string filename) {

    // Create an output filestream object
    std::ofstream myFile(filename);

    // Send data to the stream
    for (int i = 0; i < dataset.size(); ++i)
    {
        for (int j = 0; j < dataset[i].size(); ++j)
        {
            myFile << dataset[i][j];
            if (j != dataset[i].size() - 1) myFile << ","; // No comma at end of line
        }
        myFile << "\n";
    }

    // Close the file
    myFile.close();
    cout << "Successfully written data into " << filename << endl;
}


vector<vector<float>> return2DFloatVector(vector<vector<string>> vec) {
    vector< vector<float > > matFloat;

    for (size_t i = 0; i < vec.size(); ++i)
    {
        matFloat.push_back(std::vector<float>());
        for (size_t j = 0; j < vec[i].size(); j++) {
            string s = vec[i][j];
            string chars = "﻿";
            for (char c : chars) {
                s.erase(std::remove(s.begin(), s.end(), c), s.end());
            }
            if ((!s.empty()) & (s.size() > 1)) {
                matFloat[i].push_back(stof(s));
                //cout << matRows[i][j] << " ";
            }
        }
        //cout << endl;
    }

    int Nrows = matFloat.size();
    int Ncols = matFloat[0].size();
    cout << "number of rows: " << Nrows << endl;
    cout << "number of columns: " << Ncols << endl;
    return matFloat;

}

vector<vector<string>> return2DStringVector(vector<vector<float>> vec) {
    vector< vector<string > > matString;

    for (size_t i = 0; i < vec.size(); ++i)
    {
        matString.push_back(std::vector<string>());
        for (size_t j = 0; j < vec[i].size(); j++) {
            matString[i].push_back(to_string(vec[i][j]));
            //cout << matRows[i][j] << " ";
        }
        //cout << endl;
    }

    int Nrows = matString.size();
    int Ncols = matString[0].size();
    cout << "number of rows: " << Nrows << endl;
    cout << "number of columns: " << Ncols << endl;
    return matString;

}

tuple<vector<vector<float>>, vector<vector<float>>, vector<string>> returnVecwithKsize(vector<vector<float>> disVec,
    vector<vector<float>> conVec, vector<string> geneIdVec, int k) {
    
    for (size_t i = 0; i < disVec.size(); ++i)
    {
        vector<float> myVec1, myVec2;
        for (size_t j = 0; j < disVec[i].size(); j++) {
            myVec1.push_back(disVec[i][j]);
        }
        for (size_t l = 0; l < conVec[i].size(); l++) {
            myVec2.push_back(conVec[i][l]);
        }
        if ((myVec1.size() < k) | (myVec2.size() < k)) {
            disVec.erase(disVec.begin()+i);
            conVec.erase(conVec.begin()+i);
            geneIdVec.erase(geneIdVec.begin()+i);
        }
        myVec1.clear();
        myVec2.clear();
    }

    cout << "number of Diesease rows: " << disVec.size() << endl;
    cout << "number of Control rows: " << conVec.size() << endl;
    cout << "number of Gene ID rows: " << geneIdVec.size() << endl;

    return { disVec, conVec, geneIdVec };

}

// Function to find mean.
float Mean(vector<float> vec)
{
    float sum = 0;
    for (auto& elem : vec)
    {
        sum = sum + elem;
    }
    return sum / vec.size();
}

// Function to find standard
// deviation of given vector.
float standardDeviation(vector<float> vec)
{
    float sum = 0;
    for (auto& elem : vec)
    {
        sum = sum + (elem - Mean(vec)) *
            (elem - Mean(vec));
    }
    return sqrt(sum / (vec.size() - 1));
}

// Function to find t-test of
// two set of statistical data.
float tTest(vector<float> vec1, vector<float> vec2)
{
    float meanVal1 = Mean(vec1);
    float meanVal2 = Mean(vec2);
    float stdVal1 = standardDeviation(vec1);
    float stdVal2 = standardDeviation(vec2);

    // Formula to find t-test
    // of two set of data.
    float t_test = (meanVal1 - meanVal2) / sqrt((stdVal1 * stdVal2)
        / vec1.size() + (stdVal1 * stdVal2) / vec2.size());
    return t_test;
}

pair<vector<float>, vector<float>> getMeanandStdVectors(vector<vector<float>> v) {
    vector<float> meanVec;
    vector<float> stdVec;

    for (size_t i = 0; i < v.size(); ++i)
    {
        vector<float> myvector;
        for (size_t j = 0; j < v[i].size(); j++) {
            myvector.push_back(v[i][j]);
        }
        float meanVal = Mean(myvector);
        float stdVal = standardDeviation(myvector);
        meanVec.push_back(meanVal);
        stdVec.push_back(stdVal);
        myvector.clear();
    }

    return { meanVec, stdVec };
}

vector<float> getTtest(vector<vector<float>> disVec, vector<vector<float>> conVec) {
    vector<float> ttestVec;

    for (size_t i = 0; i < disVec.size(); ++i)
    {
        vector<float> diseaseVec, controlVec;
        for (size_t j = 0; j < disVec[i].size(); j++) {
            diseaseVec.push_back(disVec[i][j]);
        }
        for (size_t k = 0; k < conVec[i].size(); k++) {
            controlVec.push_back(conVec[i][k]);
        }
        float tTestVal = tTest(diseaseVec, controlVec);
        ttestVec.push_back(tTestVal);
        diseaseVec.clear();
        controlVec.clear();
    }

    return ttestVec;
}

// A function to randomly select
// k items from stream[0..n-1].
vector<float> returnKItems(vector<float> vec, int k)
{   
    int i = 0;
    random_device rd;
    mt19937 g(rd());
    srand(time(NULL));
    shuffle(vec.begin(), vec.end(), g);
    vector<float> reservoir;
    for (auto& elem : vec) {
        if (i < k) {
            reservoir.push_back(elem);
        }
        i++;
        }
    return reservoir;
}

vector<float> replaceItems(vector<float> originalVec, vector<float> sampleVec1, vector<float> sampleVec2)
{
    for (int i=0; i < sampleVec1.size(); i++)
    {
        replace(originalVec.begin(), originalVec.end(), sampleVec1[i], sampleVec2[i]);
    }
    return originalVec;
}

vector<vector<float>> randomPermutation(vector<vector<float>> disVec, vector<vector<float>> conVec) {
    vector<vector<float>> randomTtestVec;
    randomTtestVec.resize(disVec.size());

    for (size_t i = 0; i < disVec.size(); ++i)
    {
        vector<float> diseaseVec, controlVec, sampleDisease, sampleControl;
        for (size_t j = 0; j < disVec[i].size(); j++) {
            diseaseVec.push_back(disVec[i][j]);
        }
        for (size_t k = 0; k < conVec[i].size(); k++) {
            controlVec.push_back(conVec[i][k]);
        }
        for (int m = 0; m < 1000; m++) {
            srand((unsigned)time(NULL) + m);
            sampleDisease = returnKItems(diseaseVec, 3);
            sampleControl = returnKItems(controlVec, 3);
            vector<float> newDiseaseVec = replaceItems(diseaseVec, sampleDisease, sampleControl);
            vector<float> newControlVec = replaceItems(controlVec, sampleControl, sampleDisease);
            float randomtTest = tTest(newDiseaseVec, newControlVec);
            randomTtestVec[i].push_back(randomtTest);
        }
        //cout << "Index: " << i << endl;
        //cout << "Random Permutation size: " << randomTtestVec[i].size() << endl;
        diseaseVec.clear();
        controlVec.clear();
    }
    return randomTtestVec;

}

vector<float> calculateDScore(vector<vector<float>> randPermVec, vector<float> ttestVec) {
    vector<float> dscoreVec, newVec;
    vector<string> newGeneVec, sortGenes;

    for (size_t i = 0; i < randPermVec.size(); ++i)
    {
        vector<float> randPerm;
        for (size_t j = 0; j < randPermVec[i].size(); j++) {
            randPerm.push_back(randPermVec[i][j]);
        }
        float meanRandPerm = Mean(randPerm);
        float stdRandPerm = standardDeviation(randPerm);
        float dScore = (ttestVec[i] - meanRandPerm) / stdRandPerm;
        dscoreVec.push_back(dScore);
        randPerm.clear();
    }
    

    return dscoreVec;
}

pair<vector<float>, vector<string>> sortDScore(vector<float> dscoreVec, vector<string> genVec) {
    vector<float> newVec;
    vector<string> sortGenes;

    //Sort in Descending order
    newVec.assign(dscoreVec.begin(), dscoreVec.end());
    vector<size_t> newVec2 = sort_indexes(newVec);
    reverse(newVec2.begin(), newVec2.end());
    newVec.clear();

    for (auto& elem : newVec2) {
        newVec.push_back(dscoreVec[elem]);
        sortGenes.push_back(genVec[elem]);
        cout << "D-Score: " << dscoreVec[elem] << "\tGene Index: " << genVec[elem] << endl;
        
    }
    return { newVec, sortGenes };
}


int main()
{
    vector< vector<float > > matControl, matRenal, randomPermVecs;
    vector< vector<string > > controlTwodVec, renalTwodVec, geneTwodVec, randPermStrVec;
    vector<string> geneVec, sortGeneVec;
    //matControl = parseCsvToVector("../control.csv");
    renalTwodVec = parseCsvToVector("renal.csv");
    controlTwodVec = parseCsvToVector("others.csv");
    geneTwodVec = parseCsvToVector("gene_index.csv");
    geneVec = flatten(geneTwodVec);
    matRenal = return2DFloatVector(renalTwodVec);
    matControl = return2DFloatVector(controlTwodVec);
    auto tuples = returnVecwithKsize(matRenal, matControl, geneVec, 3); //remove gene with less than 3 patients
    matRenal = get<0>(tuples);
    matControl = get<1>(tuples);
    geneVec = get<2>(tuples);
    vector<float> meanRenal, stdRenal, meanControl, stdControl, ttstVec, dScoreVec;
    auto renalVecs = getMeanandStdVectors(matRenal);
    meanRenal = renalVecs.first;
    stdRenal = renalVecs.second;
    auto conVecs = getMeanandStdVectors(matControl);
    meanControl = conVecs.first;
    stdControl = conVecs.second;
    for (size_t i = 0; i < meanRenal.size(); ++i)
    {
        cout << "Gene Index: " << geneVec[i] << endl;
        cout << "mean Renal: " << meanRenal[i] << "\t mean Control: " << meanControl[i] << endl;
        cout << "Std Renal: " << stdRenal[i] << "\t std Control: " << stdControl[i] << endl;
    }
    auto start = std::chrono::high_resolution_clock::now();
    ttstVec = getTtest(matRenal, matControl);
    /*for (size_t i = 0; i < ttstVec.size(); ++i)
    {
        cout << "Gene Index: " << geneVec[i] << endl;
        cout << "T test value: " << ttstVec[i] << endl;
    }*/
    randomPermVecs = randomPermutation(matRenal, matControl);
    //randPermStrVec = return2DStringVector(randomPermVecs);
    //write_csv(randPermStrVec, "random_perm.csv");
    dScoreVec = calculateDScore(randomPermVecs, ttstVec);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop -start);
    printf("\nTotal Execution time %ld seconds\n", duration.count());
    auto sortVecs = sortDScore(dScoreVec, geneVec);
    dScoreVec = sortVecs.first; 
    sortGeneVec = sortVecs.second;
    geneTwodVec.clear();
    geneTwodVec.resize(sortGeneVec.size());
    for (size_t i = 0; i < geneTwodVec.size(); ++i) {
        geneTwodVec[i].push_back(sortGeneVec[i]);
        geneTwodVec[i].push_back(to_string(dScoreVec[i]));
    }
    write_csv(geneTwodVec, "sort_gene_dscore.csv");

    return 0;

}