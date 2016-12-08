#include<iostream>
#include "sas.h"

#include <iomanip>
#include <random>
#include <cmath>
#include <sys/time.h>
#include <ctime>

const double PI = 3.14159265359;

std::random_device rd;
std::mt19937 gen(rd());

std::normal_distribution<double> ___normal_distribution(0.0, 1.0);

double _normal_distribution(double mean, double dev, double x){
    
    return 1.0/(sqrt(2 * PI) * dev) * exp((-1.0/2.0) * pow((x - mean) / dev, 2));
}


/* generate a batch of random number at the same time */
void __normal_distribution(int samplecount, std::vector<double> &samples){
    for(int i = 0; i <samplecount; i++) {
        samples[i] = ___normal_distribution(gen);
    }
}

double cramer_von_mises_iid(int bincount, int samplecount){
    
    std::vector<double> pmf;
    std::vector<double> values;
    
    double h = (20.0) / (bincount - 1); //[take value in [-4, 4]
    for(int i = 0; i < bincount; i++){
        double t = -10 + i * h;
        values.push_back(t);
        pmf.push_back(_normal_distribution(0.0, 1.0, t) + 0.02);
    }
    
    SystematicAliasSampling _sample(pmf, values);
    _sample.aliastable();

    std::vector<double> samples(samplecount); // store the results
    std::vector<double> cdf(bincount);
    _sample.cummulative_distribution(cdf);
    
    std::uniform_real_distribution<> __random(0.0, 1.01);
    
    /* binary search */
    for(int i = 0; i < samplecount; i++){
        double x  = __random(gen);
        if(x > 1.0) x = 1.0;
        size_t index = std::upper_bound(cdf.begin(), cdf.end(), x) - cdf.begin();
        samples[i] = values[index];
    }
    
    
    std::vector<double> edf(bincount);
   
    _sample.empirical_distribution(samples, edf);
   
    double cvm = _sample.cramer_von_mises(cdf, edf);
    
    return cvm;
}
double cramer_von_mises_sas(int bincount, int samplecount){
    
    std::vector<double> pmf;
    std::vector<double> values;
    
    double h = (20.0) / (bincount - 1); //[take value in [-4, 4]
    for(int i = 0; i < bincount; i++){
        double t = -10 + i * h;
        values.push_back(t);
        pmf.push_back(_normal_distribution(0.0, 1.0, t) +0.02);
    }
    
    SystematicAliasSampling _sample(pmf, values);
    
    std::vector<double> samples(samplecount); // store the results
    _sample.aliastable(); //generate table
    
    std::vector<double> cdf(bincount);
    _sample.cummulative_distribution(cdf);
    

    //_sample.systematicaliassampling(samplecount, samples);
    _sample.simple_sas(samplecount, samples);
    //_sample.goldratioaliassampling(samplecount, samples);

    std::vector<double> edf(bincount);
    _sample.empirical_distribution(samples, edf);
    
    double cvm = _sample.cramer_von_mises(cdf, edf);
    
    return cvm;
}
int main(int argc, char *argv[]){

    
    int samplecount = 0;
    int bincount = 0;

    if(argc != 3){
        std::cout << "Usage: ./test bincount samplecount" << std::endl;
        exit(1);
    } 
    bincount = atoi(argv[1]); //bincount
   // samplecount = atoi(argv[2]); //sampple count
    
    /*
    struct timeval tvs, tvm; 
    std::vector<double> pmf;
    std::vector<double> values;
    
    double h = (8.0) / (bincount - 1); //[take value in [-4, 4]
    
    for(int i = 0; i < bincount; i++){
        double t = -4 + i * h;
        values.push_back(t);
        pmf.push_back(_normal_distribution(0.0, 1.0, t));
    }
    
    //  initial sas class 
    SystematicAliasSampling _sample(pmf, values);
    
    
    std::vector<double> samples(samplecount); // store the results
    _sample.aliastable(); //generate table
    //std::vector<double> av = _sample.aliasvalue;;
    //std::sort(av.begin(), av.end());

    //for(int i = 0 ; i < av.size(); i++){
        //std::cout << av[i] << std::endl;
   //     std::cout << _normal_distribution(0.0, 1.0, av[i]) << std::endl;
    //} 
    std::vector<double> cdf(bincount);
    _sample.cummulative_distribution(cdf);
    
    std::uniform_real_distribution<> __random(0.0, 1.01);
    
    gettimeofday(&tvs, NULL);
    
    // binary search 
    std::vector<double> b_samples(samplecount);
    for(int i = 0; i < samplecount; i++){
        double x  = __random(gen);
        if(x > 1.0) x = 1.0;
        b_samples[i] = values[std::upper_bound(cdf.begin(), cdf.end(), x) - cdf.begin()];
        //std::cout << std::upper_bound(cdf.begin(), cdf.end(), x) - cdf.begin() << std::endl;
        
    }
    gettimeofday(&tvm, NULL);
    double span0 = tvm.tv_sec-tvs.tv_sec + (tvm.tv_usec-tvs.tv_usec)/1000000.0;
    std::cout << "Binary Serach cost: " << span0 <<std::endl;
    
    // alias sampling method 
    std::uniform_real_distribution<double> ___random(0.0, bincount);
    gettimeofday(&tvs, NULL);
    std::vector<double> alias_samples(samplecount);
    for(int i = 0; i < samplecount; i++){
        double x  = ___random(gen);
        alias_samples[i]  = _sample.aliassample(int(x), x - int(x));
    }
    gettimeofday(&tvm, NULL);
    double span = tvm.tv_sec-tvs.tv_sec + (tvm.tv_usec-tvs.tv_usec)/1000000.0;
    std::cout << "Alias method cost: " << span <<std::endl;
   

    // sas sampling 
    gettimeofday(&tvs, NULL);
    _sample.systematicaliassampling(samplecount, samples);
    gettimeofday(&tvm, NULL);
    
    
    double span1 = tvm.tv_sec-tvs.tv_sec + (tvm.tv_usec-tvs.tv_usec)/1000000.0;
    std::cout << "SAS cost: " << span1 <<std::endl;
    
    //sas golden sampling 
    gettimeofday(&tvs, NULL);
    std::vector<double> samples2(samplecount);
    _sample.goldratioaliassampling(samplecount, samples2);
    
    gettimeofday(&tvm, NULL);
    double span4 = tvm.tv_sec-tvs.tv_sec + (tvm.tv_usec-tvs.tv_usec)/1000000.0;
    std::cout << "SAS-golden cost: " << span4 <<std::endl;
   
    // normal_distribution in c++11 
    gettimeofday(&tvs, NULL);
    std::vector<double> samples3(samplecount);
    __normal_distribution(samplecount, samples3);
    gettimeofday(&tvm, NULL);
    double span2 = tvm.tv_sec-tvs.tv_sec + (tvm.tv_usec-tvs.tv_usec)/1000000.0;
    std::cout << "normal_distribution cost: " << span2 <<std::endl;
    
    //discrete_distribution sampling in c++11 
    
    std::discrete_distribution<> _discrete_distribution(pmf.begin(), pmf.end());
    gettimeofday(&tvs, NULL);
    std::vector<double> __p = _discrete_distribution.probabilities();
    for(int i = 0 ; i < samplecount; i++){
        _discrete_distribution(gen);
    }
    gettimeofday(&tvm, NULL);
    double span3 = tvm.tv_sec-tvs.tv_sec + (tvm.tv_usec-tvs.tv_usec)/1000000.0;
    std::cout << "discrete_distribution cost: " << span3 <<std::endl;
    

   // std::vector<double> edf(bincount);
   
    //_sample.empirical_distribution(b_samples, edf);
    
    */
   // for(int i = 0; i < bincount ; i++){
    //    std::cout << edf[i]<< std::endl;
        
   // }
    for(int j = 0 ; j < 1000; j ++){
        for(int i = 20; i < 2 * bincount; i++){
           // std::cout << "+++++++++++++"; 
           // double cvm = cramer_von_mises_alias(bincount, i);
            double cvm = cramer_von_mises_sas(bincount, i);
            //double cvm = cramer_von_mises_iid(bincount, i);
            std::cout << cvm <<" ";
        }
        std::cout << std::endl;
    }
    return 0; 

}

