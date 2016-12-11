#ifndef _SAS_H_ 
#define _SAS_H_

/*This the implementation of the systematic alias sampling method,
 * which is the combination of the alias method and systematic method. 
 *
 * This method can generate k samplings at the same in O(k) time cost 
 *
 * The main references for this program is inspired by the paper: Systematic Alias
  Sampling: An Efficient and Low-Variance Way  to Sample from a Discrete Distribution written by I  LARI VALLIVAARA, KATJA POIKSELKA , PAULI RIKULA, and JUHA RO NING 
 * 
*/

#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>
#include <stdio.h>
#include <string.h>
#include <string>
#include <cmath>
#include <random>
#include <stack>
#include <algorithm>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> __random(0.0, 1.0);

class SystematicAliasSampling{

    private:

        typedef std::vector<double> Table; //Typedef the vector<double> type  
        typedef std::vector<double> Values; 
        /* initial value table */
        Table pmf; // probability mass function;
        Values values;  // values corresponding to the probability mass function
        /* initial value table */
    //private:
        /* Tables after constructing the Alias Table */
        Table aliasvalue;
        Table aliasindices;
        Table aliasprobabilities;
        
        size_t bincount;
        /* Tables after constructing the Alias Table */
        
        double eps = 0.07; // parameter to determine the divisiable;

        int minbatchsize = 16;
        int minrecursize = 4 * minbatchsize;
        //int minrecursize = 4 * minbatchsize;
        
        int batchsplitnumerator = 7;
        int batchsplitdenominator = 13;
        

    public:
        /* empty constructor */
        SystematicAliasSampling(){} 
        

        /* override  constructor function*/
        SystematicAliasSampling(const Table & _pmf, const Values &_value){
            bincount = _pmf.size();
            pmf = _pmf;
            values = _value;
        }

        /* initial function to initize the value */
        void init(const Table & _pmf, const Values &_value){
            bincount = _pmf.size();
            pmf = _pmf;
            values = _value;
        }



    private:
        
        /* normalize the arrary */
        bool normalize(){ 
            double sum = 0.0;
            for(auto p: pmf) sum += p; //get the sum for array
            for(auto &p: pmf) p = p / sum; //normalize the array
            sum = 0.0;
            for(auto p:pmf) sum += p;
            return (fabs(sum - 1.0) < 1e-6); //make sure the sum is close to 1.0
        }
    public:
        void aliastable(){
            
            if( !normalize() ) std::cout << "Wrong normalization" <<std::endl;;
            Table F(bincount); //define F
            for(size_t i = 0; i < bincount ; i++){
                F[i] = pmf[i] * double(bincount);
            }
            aliasvalue.resize(bincount);
            aliasindices.resize(bincount);
            aliasprobabilities.resize(bincount);
            std::stack<double> G, S;

            /* Denote by G the set of i, if F[i] >= 1, otherwise S*/
            
            for(size_t i = 0; i < bincount ; i++) {
                if(F[i] >= 1.0 || fabs(F[i] - 1.0) < 1e-6) G.push(i);
                else S.push(i);
                
                aliasindices[i] = i;
            }
             
            /* main loop until S is empty */
            while(!G.empty() && !S.empty()){
            
               int j = S.top();
               S.pop(); // remove the top element of stack S
               int k = G.top(); //get the top element in the stack G
               aliasindices[j] = k;
               aliasprobabilities[j] = 1.0 - F[j];
               
               F[k] = (F[k] + F[j]) -1.0;
               if(F[k] < 1.0){
                    G.pop();
                    S.push(k);
               }

            }
            for(size_t i = 0; i < bincount ; i++) {
                aliasvalue[i] = values[aliasindices[i]];
            }
         //   aliasprobabilities = F;
        }
        
    private:    
        bool aldivisiable(size_t _bincount, size_t ssize){

            double q = _bincount / double(ssize);
            double frac = q - int(q);
            return std::min(frac, 1.0 - frac) < eps;
        }
        
        bool isdivisible(size_t _bincount, size_t ssize){
            return (aldivisiable(_bincount, ssize) ||aldivisiable(_bincount *4, ssize) || 
                aldivisiable(_bincount * 5, ssize) || aldivisiable(_bincount * 6, ssize));
        }
    public:

        /* alias method: generate one sampling each time */
        double aliassample(int intpart, double fracpart){
            
            if(fracpart <= aliasprobabilities[intpart]){
                return aliasvalue[intpart];
                
                //return values[intpart];
            }else{
                return values[intpart];
                //return aliasvalue[intpart];
            }
        }

        /* systematicaliassampling implemetation: Using samples to store the
         * results */ 
        void systematicaliassampling(int samplecount, Table &samples, int from = 0){
            int splitindex;  
            if(samplecount > minbatchsize && isdivisible(bincount, samplecount)){
                if(samplecount <= minrecursize){
                    splitindex = minbatchsize;
                }else{
                    splitindex = samplecount * batchsplitnumerator / batchsplitdenominator;
                }
                systematicaliassampling(samplecount - splitindex, samples, fillfrom);
                systematicaliassampling(splitindex, samples, fillfrom + samplecount - splitindex);
            }else{
                 
                double steps = double(bincount) / samplecount;
                double r = steps * (1.0 - __random(gen));
                //double x =  r;
                double x =  bincount - r;
                int i = from;

                int to = from + samplecount;
                while(i < to){
                    int ri = int(x);
                    double rf  = x  -ri;
                    samples[i] = aliassample(ri, rf);
                    //x += steps;
                    x -= steps;
                    i++;
                }
            }
        
        }

        /* sas method without taking care of divisible problem */
        void simple_sas(int samplecount, Table &samples, int fillfrom = 0){
                
                double steps = double(bincount) / samplecount;
                //double r = _random(0.0, steps);
                double r = steps * (1.0 - __random(gen));
                double x =  bincount - r;
                int i = fillfrom;

                int fillto = fillfrom + samplecount;
                while(i < fillto){
                  //  int ri = int(x);
                   // int rf  = x  -ri;
                    samples[i] = aliassample(int(x), x- int(x));
                    x -= steps;
                    i++;
                }
        }
        /*Gold ratio alias sampling: results in samples */
        void goldratioaliassampling(int samplecount, Table &samples){
            double x = __random(gen);
            double ratio = 0.618033988749895;
            int i = 0;
            while(i < samplecount){
                x += ratio;
                x -= int(x);
                double d = bincount * x;
                samples[i] = aliassample(int(d), d - int(d));
                i++;
            }
        }

        void empirical_distribution(Table &samples, Values &edf){
            std::sort(samples.begin(), samples.end()); 
            size_t count = 0;
            double last = samples[samples.size() - 1];
            for(size_t i = 0; i < values.size(); i++){
                count = 0;
                if(values[i] < samples[0] && fabs(samples[0] - values[i]) > 1e-6){
                    edf[i] = 0.0;
                    continue;
                }else if(values[i] >= last || fabs(values[i] - last)  < 1e-6){
                
                    edf[i] = 1.0;
                    continue;
                }else{
                    for(size_t j = 0 ; j < samples.size(); j++){
                        if( values[i] >= samples[j] || fabs(values[i] - samples[j]) < 1e-6 ){
                            count++;
                        }else{
                            break;
                        }

                    }
                    edf[i] = double(count) / samples.size();
                }
            }
        
        }

        void cummulative_distribution(Table &cdf){
            for(size_t i = 0; i < bincount; i++){
                if(i == 0) cdf[i] = pmf[i];
                else{
                    cdf[i] = cdf[i-1] + pmf[i];
                }
            }
        }

        double cramer_von_mises(const Table &cdf, const Table &edf){
            double sum = 0.0;
            for(size_t i = 0; i < bincount; i++){
                sum += (edf[i] - cdf[i]) * (edf[i] - cdf[i]);
            }
            return sqrt(sum / bincount);
        }



};

#endif
