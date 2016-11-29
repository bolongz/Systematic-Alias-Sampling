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


class SystematicAliasSampling{

    public:

        typedef std::vector<double> Table; //Typedef the vector<double> type  
        typedef std::vector<double> Values; 
        /* initial value table */
        Table pmf; // probability mass function;
        Values values;  // values corresponding to the probability mass function
        /* initial value table */
    private:
        /* Tables after constructing the Alias Table */
        Table aliasvalue;
        Table aliasindices;
        Table aliasprobabilities;
        
        size_t bincount;
        /* Tables after constructing the Alias Table */
        
        double eps = 0.07; // parameter to determine the divisiable;

        int minbatchsize = 16;
        int minrecursize = 4 * minbatchsize;
        
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
        bool normalize(Table &_pmf){ 
            double sum = 0.0;
            for(auto p: _pmf) sum += p; //get the sum for array
            for(auto &p: _pmf) p = p / sum; //normalize the array
            sum = 0.0;
            for(auto p:_pmf) sum += p;
            return (fabs(sum - 1.0) < 1e-6); //make sure the sum is close to 1.0
        }
    public:
        void aliastable(){
            
            normalize(pmf);
            //size_t n = pmf.size();
            Table F(bincount); //define F
            for(size_t i = 0; i < bincount ; i++){
                F[i] = pmf[i] * bincount;
                //std::cout << pmf[i] << " " << F[i] << std::endl;
            }
            aliasvalue.resize(bincount);
            aliasindices.resize(bincount);
            aliasprobabilities.resize(bincount);
            std::stack<double> G, S;

            /* Denote by G the set of i, if F[i] >= 1, otherwise S*/
            
            for(size_t i = 0; i < bincount ; i++) {
                if(F[i] >= 1.0) G.push(i);
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
        }
    public:
        double _random(double _min, double _max){
            std::random_device rd;

            /* generate the random double value from [0.0, bincount) */
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> __random(_min, _max);
            return __random(gen);
        }
        
    private:    
        bool aldivisiable(size_t _bincount, size_t ssize){

            double q = _bincount / double(ssize);
            double frac = q - int(q);
            return std::min(frac, 1.0 - frac) < eps;
        }
        
        bool isdivisible(size_t _bincount, size_t ssize){
            return aldivisiable(_bincount, ssize) ||aldivisiable(_bincount *4, ssize) || 
                aldivisiable(_bincount * 5, ssize) || aldivisiable(_bincount * 6, ssize);
        }
    public:

        /* alias method: generate one sampling each time */
        double aliassample(int intpart, double fracpart){
            
            if(fracpart <= aliasprobabilities[intpart]){
                return aliasvalue[intpart];
            }else{
                return values[intpart];
            }
        }

        /* systematicaliassampling implemetation: Using samples to store the
         * results */ 
        void systematicaliassampling(int samplecount, Table &samples, int fillfrom = 0){
            int splitindex;  
            if(samplecount > minbatchsize && isdivisible(bincount, samplecount)){
                if(samplecount <= minrecursize){
                    splitindex = minbatchsize;
                }else{
                    splitindex = int(samplecount * (double(batchsplitnumerator) 
                                / batchsplitdenominator));
                
                }
                systematicaliassampling(samplecount - splitindex, samples, fillfrom);
                systematicaliassampling(splitindex, samples, fillfrom + samplecount - splitindex);
            }else{
                
                double steps = double(bincount) / samplecount;
                double r = steps * (1.0 - _random(0.0, 1.0));
                double x =  bincount - r;
                int i = fillfrom;

                int fillto = fillfrom + samplecount;
                while(i < fillto){
                    int ri = int(x);
                    int rf  = x  -ri;
                    samples[i] = aliassample(ri, rf);
                    x -= steps;
                    i++;
                }
            }
        
        }
        /*Gold ratio alias sampling: results in samples */
        void goldratioaliassampling(int samplecount, Table &samples){
            double x = _random(0.0, 1.0);
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



};

#endif
