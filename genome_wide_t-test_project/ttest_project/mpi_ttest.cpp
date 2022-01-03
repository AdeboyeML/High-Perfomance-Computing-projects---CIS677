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
#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
#define SIZE 4532
#define renalSIZE 8
#define controlSIZE 52
#define randSize 5000

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

void return_2d_renal_array(vector<vector<float>> vec, float arr[SIZE][renalSIZE]) {

    //cout << "2D Float Arrays: " << endl;
    for (size_t i = 0; i < vec.size(); ++i)
        {
            //arr[i] = new float[vec[i].size()];
            for (size_t j = 0; j < vec[i].size(); j++) {
                arr[i][j] = vec[i][j];
                //cout << arr[i][j] << " ";
            }
            //arrFloat[i] = new float[vec[i].size()];
            
            //cout << endl;
        }

    /*int Nrows = sizeof (arr) / sizeof (arr[0]);
    int Ncols = sizeof (arr[0]) / sizeof (int);
    cout << "number of rows: " << Nrows << endl;
    cout << "number of columns: " << Ncols << endl;*/
}

void return_2d_ctrl_array(vector<vector<float>> vec, float arr[SIZE][controlSIZE]) {

    //cout << "2D Float Arrays: " << endl;
    for (size_t i = 0; i < vec.size(); ++i)
        {
            //arr[i] = new float[vec[i].size()];
            for (size_t j = 0; j < vec[i].size(); j++) {
                arr[i][j] = vec[i][j];
                //cout << arr[i][j] << " ";
            }
            //arrFloat[i] = new float[vec[i].size()];
            
            //cout << endl;
        }

    /*int Nrows = sizeof (arr) / sizeof (arr[0]);
    int Ncols = sizeof (arr[0]) / sizeof (int);
    cout << "number of rows: " << Nrows << endl;
    cout << "number of columns: " << Ncols << endl;*/
}


vector<float> return1DFloatVec(float arr[SIZE][1]) {
    vector<float> matFloat;

    for (size_t i = 0; i < SIZE; ++i)
    {
        matFloat.push_back(arr[i][0]);
    }

    int Nrows = matFloat.size();
    cout << "number of rows: " << Nrows << endl;
    return matFloat;

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
float Mean(float arr[], int n)
{
    float sum = 0;
    for (int i = 0; i < n; i++)
        sum = sum + arr[i];
    return sum / n;
}
 
// Function to find standard
// deviation of given array.
float standardDeviation(float arr[], int n)
{
    float sum = 0;
    for (int i = 0; i < n; i++)
        sum = sum + (arr[i] - Mean(arr, n)) *
                    (arr[i] - Mean(arr, n));
 
    return sqrt(sum / (n - 1));
}
 
// Function to find t-test of
// two set of statistical data.
float tTest(float arr1[], int n,
            float arr2[], int m)
{
    float mean1 = Mean(arr1, n);
    float mean2 = Mean(arr2, m);
    float sd1 = standardDeviation(arr1, n);
    float sd2 = standardDeviation(arr2, m);
 
    // Formula to find t-test
    // of two set of data.
    float t_test = (mean1 - mean2) / sqrt((sd1 * sd1)
                              / n + (sd2 * sd2) / m);
    return t_test;
}


float *replaceItems(float *arr, float *oldArr, float *newArr)
{
    for (int i=0; i < sizeof (oldArr) / sizeof (oldArr[0]); i++)
    {
        replace(arr, arr + (sizeof (arr) / sizeof (arr[0])), oldArr[i], newArr[i]);
    }
    return arr;
}


int main(int argc, char **argv)
{
    int from, to;
    float renal2DArray [SIZE][renalSIZE]; 
    float control2DArray [SIZE][controlSIZE];
    float dScoreVec [SIZE][1];
    fill(renal2DArray[0], renal2DArray[0] + SIZE * renalSIZE, 0.000);
    fill(control2DArray[0], control2DArray[0] + SIZE * controlSIZE, 0.000);
    fill(dScoreVec[0], dScoreVec[0] + SIZE * 1, 0.000);
    vector< vector<float > > matControl, matRenal;
    vector< vector<string > > controlTwodVec, renalTwodVec, geneTwodVec, randPermStrVec;
    vector<string> geneVec, sortGeneVec;
    vector<float>ttstVec, dScores, newVec;
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
    return_2d_renal_array(matRenal, renal2DArray);
    return_2d_ctrl_array(matControl, control2DArray);

    /*for (size_t i = 0; i < SIZE; ++i)
    {
        //Ncols = matControl[i].size();
        cout << i << " <--Index NUmber"<<  endl;
        int m = 0;
        //cout << Ncols << " <--num of columns"<<  endl;
        for (size_t j = 0; j < controlSIZE; j++) {
            if (abs(control2DArray[i][j]) > 0.000) {
                cout << control2DArray[i][j] << " ";
                m++;
            }
            //cout << control2DArray[i][j] << " ";
        }
        

        cout << "\n Column Size: " << m << endl;
    }*/
    
    MPI_Init(NULL, NULL);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    from = world_rank * SIZE/world_size;
    to = ((world_rank+1) * SIZE/world_size);

    /*if(world_rank==0){    
        
    }*/
    double start_time = MPI_Wtime();
    MPI_Bcast (renal2DArray, SIZE*renalSIZE, MPI_INT, 0, MPI_COMM_WORLD);
    //printf("\n\n");
    if(world_rank==0){
        MPI_Scatter (&control2DArray[0][0], SIZE*controlSIZE/world_size, MPI_INT, MPI_IN_PLACE, SIZE*controlSIZE/world_size, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else{
        MPI_Scatter (&control2DArray[0][0], SIZE*controlSIZE/world_size, MPI_INT, &control2DArray[from][0], SIZE*controlSIZE/world_size, MPI_INT, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("computing processors %d (from row %d to %d)\n",world_rank, from, to-1);
    for (size_t i=from; i<to; i++) {
        int ren = 0;
        float subRenal[renalSIZE];
        for (size_t k=0; k<renalSIZE; k++){
            if (abs(renal2DArray[i][k]) > 0.000) {
                subRenal[ren] = renal2DArray[i][k];
                ren++;
            }
        }
        int ctr = 0;
        float subControl[controlSIZE];
        for (size_t v=0; v<controlSIZE; v++){
            if (abs(control2DArray[i][v]) > 0.000) {
                subControl[ctr] = control2DArray[i][v];
                ctr++;
            }
        }
        float tTst= tTest(subRenal, ren, subControl, ctr);
        float randomPerms[randSize];
        for (int n = 0; n < randSize; n++) {
            srand((unsigned)time(NULL) + world_rank*world_size + n);
            random_device rd;
            mt19937 g(rd());
            shuffle(subRenal, subRenal + ren, g);
            shuffle(subControl, subControl + ctr, g);
            float renalRand[3] = {subRenal[0], subRenal[1], subRenal[2]};
            float ctrlRand[3] = {subControl[0], subControl[1], subControl[2]};
            float *newRenal = replaceItems(subRenal, renalRand, ctrlRand);
            float *newCtrl = replaceItems(subControl, ctrlRand, renalRand);
            randomPerms[n] = tTest(newRenal, ren, newCtrl, ctr);
        }
        dScoreVec[i][0] = (tTst - Mean(randomPerms, randSize))/standardDeviation(randomPerms, randSize);
        fill(subRenal, subRenal+renalSIZE, 0);
        fill(subControl, subControl+controlSIZE, 0);
        fill(randomPerms, randomPerms+randSize, 0);
        
    }
    
    if(world_rank==0){
         MPI_Gather (MPI_IN_PLACE, SIZE/world_size, MPI_INT, &dScoreVec[0][0], SIZE/world_size, MPI_INT, 0, MPI_COMM_WORLD);
    }else{
         MPI_Gather (&dScoreVec[from][0], SIZE/world_size, MPI_INT, &dScoreVec[0][0], SIZE/world_size, MPI_INT, 0, MPI_COMM_WORLD);
    }
    double end_time = MPI_Wtime();

    if (world_rank==0) {
        //sort in descending order
        dScores = return1DFloatVec(dScoreVec);
        newVec.assign(dScores.begin(), dScores.end());
        vector<size_t> newVec2 = sort_indexes(newVec);
        reverse(newVec2.begin(), newVec2.end());
        newVec.clear();
        geneTwodVec.clear();
        geneTwodVec.resize(SIZE);
        printf("\n\n");
        int i = 0;
        {
            for (auto& elem : newVec2) {
                cout << "Gene Index: " << geneVec[elem] << endl;
                printf("D score: %f ", dScores[elem]);
                printf("\n");
                geneTwodVec[i].push_back(geneVec[elem]);
                geneTwodVec[i].push_back(to_string(dScores[elem]));
                i++;
            }
            printf("Total Execution time: %f seconds \n", end_time - start_time);
            write_csv(geneTwodVec, "sort_gene_dscore_MPI.csv");
        }

        printf("\n\n");
       
    }

    MPI_Finalize();
    return 0;

}