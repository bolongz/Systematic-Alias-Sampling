#include<iostream>
#include "sas.h"

#include <iomanip>
#include <random>
#include <cmath>

const double PI = 3.14159265359;


double _normal_distribution(double mean, double dev, double x){
    
    return 1.0/(sqrt(2 * PI) * dev) * exp((-1.0/2.0) * pow((x - mean) / dev, 2));
}

int main(int argc, char *argv[]){

    int samplecount = 0;
    int bincount = 0;

   // if(argc == 1){
   //     samplecount = 16; //default sample count is 16
   // }else{
    if(argc != 3){
        std::cout << "Usage: ./test bincount samplecount" << std::endl;
        exit(1);
    } 
    bincount = atoi(argv[1]); //bincount
    samplecount = atoi(argv[2]); //sampple count
    //}

    //std::random_device rd;
    //std::mt19937 gen(rd());
    //std::normal_distribution<double> _normal_distribution(0.0, 1.0);
    
    std::vector<double> pmf;
    std::vector<double> values;
    double h = (8.0) / bincount; //[take value in [-4, 4]
    for(int i = 0; i < bincount + 1; i++){
        double t = -4 + i * h;
        values.push_back(t);
        pmf.push_back(_normal_distribution(0.0, 1.0, t));
    }

    SystematicAliasSampling _sample(pmf, values);
    
    _sample.aliastable(); //generate table

    
    std::vector<double> samples(samplecount);
    _sample.systematicaliassampling(samplecount, samples);
    //std::vector<double> samples2(samplecount);
   
    //_sample.goldratioaliassampling(samplecount, samples2);
    
    for(int i = 0; i < samplecount ; i++){
        std::cout << samples[i]<< std::endl;
    }

    return 0; 

}

