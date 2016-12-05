#include<iostream>
#include "sas.h"

#include <iomanip>
#include <random>
#include <cmath>
#include <sys/time.h>
#include <ctime>

const double PI = 3.14159265359;


double _normal_distribution(double mean, double dev, double x){
    
    return 1.0/(sqrt(2 * PI) * dev) * exp((-1.0/2.0) * pow((x - mean) / dev, 2));
}
    
void __normal_distribution(int samplecount){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> _normal_distribution(0.0, 1.0);
    for(int i = 0; i <samplecount; i++) {
        _normal_distribution(gen);
    }
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
    struct timeval tvs, tvm; 
    std::vector<double> pmf;
    std::vector<double> values;
    double h = (8.0) / bincount; //[take value in [-4, 4]
    for(int i = 0; i < bincount + 1; i++){
        double t = -4 + i * h;
        values.push_back(t);
        pmf.push_back(_normal_distribution(0.0, 1.0, t));
    }
    SystematicAliasSampling _sample(pmf, values);
    
    //for(int i = 0; i < bincount +1 ; i++){
    
    //    std::cout << _sample.values[i] << std::endl;
    
    //}
    
    std::vector<double> samples(samplecount);
    
    _sample.aliastable(); //generate table
    gettimeofday(&tvs, NULL);
    _sample.systematicaliassampling(samplecount, samples);
    gettimeofday(&tvm, NULL);
    
    double span1 = tvm.tv_sec-tvs.tv_sec + (tvm.tv_usec-tvs.tv_usec)/1000000.0;
    std::cout << "SAS cost: " << span1 <<std::endl;
    
    gettimeofday(&tvs, NULL);
    __normal_distribution(samplecount);
    gettimeofday(&tvm, NULL);
    double span2 = tvm.tv_sec-tvs.tv_sec + (tvm.tv_usec-tvs.tv_usec)/1000000.0;
    std::cout << "normal_distribution cost: " << span2 <<std::endl;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> _discrete_distribution(pmf.begin(), pmf.end());
    gettimeofday(&tvs, NULL);
   
    std::vector<double> __p = _discrete_distribution.probabilities();

    /*for(int i = 0 ; i <bincount + 1; i++){
    
        std::cout << pmf[i] << " "<< _sample.pmf[i] << " "<< __p[i] << std::endl; 
    }
    */
    for(int i = 0 ; i < samplecount; i++){
        _discrete_distribution(gen);
        
    }
    gettimeofday(&tvm, NULL);
    double span3 = tvm.tv_sec-tvs.tv_sec + (tvm.tv_usec-tvs.tv_usec)/1000000.0;
    std::cout << "discrete_distribution cost: " << span3 <<std::endl;
    //std::vector<double> samples2(samplecount);
    //_sample.goldratioaliassampling(samplecount, samples2);
    /*
    std::vector<double> edf(bincount + 1);
    std::vector<double> cdf(bincount + 1);

    _sample.empirical_distribution(samples, edf);
    _sample.cummulative_distribution(cdf);
    std::sort(samples2.begin(), samples2.end());
    */
    for(int i = 0; i < samplecount ; i++){
      //  std::cout << samples[i]<< std::endl;
    }
    // for(int i = 0; i < bincount + 1 ;i++){
    //    std::cout << cdf[i] << std::endl;
    //}
    //
    //double cvm = _sample.cramer_von_mises(cdf, edf);
    //std::cout << cvm <<std::endl;
    return 0; 

}

