//
//  main.cpp
//  COMP528 Assessment 1
//
//  Created by Elliot Adderton on 10/03/2014.
//

// I've commented most of the document but it should be self explanatory

#include "main.h"

// Size of the arrays
const size_t vec_size = 500000;

int main() {
    
	//memory allocation for the vectors
	double *arrayA = (double*) malloc(vec_size*sizeof(double));
	double *arrayB = (double*) malloc(vec_size*sizeof(double));
    
	//initialize all the vectors
	for(size_t i = 0; i < vec_size; i++) {
		arrayA[i] = sin(i);
		arrayB[i] = cos(i);
	}
    
    // Start the counter for serial
    tick_count serialStart = tick_count::now();
    
    // Calculate the serial value
    double serialValue = pearsonCorrelationCoefficient(arrayA, arrayB);
    
    // End the counter for serial
    tick_count serialFinish = tick_count::now();
    
    // Using fixed values when printing out
    cout << fixed;
    cout << "Serial Time: " << (serialFinish-serialStart).seconds() << endl;
    cout << "Serial Value: " << serialValue << endl;
    
    // Automatically set up the task scheduler
    task_scheduler_init init;
    // Manuallt set uo the task scheduler
    //task_scheduler_init init(4);
    
    // Start the counter for parallel
    tick_count parallelStart = tick_count::now();
    
    // Calculate the parallel value
    double parallelValue = parallelPearsonCorrelationCoefficient(arrayA, arrayB);
    
    // End the counter for parallel
    tick_count parallelFinish = tick_count::now();
    
    // House keeping !!! VERY IMPORTANT !!!
    free(arrayA);
    free(arrayB);
    
    // Print the findings
    cout << "Parallel Time: " << (parallelFinish-parallelStart).seconds() << endl;
    cout << "Parallel Value: " << parallelValue << endl;
    
    cout << "Speedup from parallelism is " << (serialFinish - serialStart).seconds()/(parallelFinish - parallelStart).seconds() << endl;
    
    // Latex table output
    //cout << vec_size << " & " << (serialFinish-serialStart).seconds() << " & " << (parallelFinish-parallelStart).seconds() << " & " << (serialFinish - serialStart).seconds()/(parallelFinish - parallelStart).seconds() << "\\\\" << endl;
}

double parallelMean(double *a) {
    return parallel_reduce(blocked_range<double*>(a,a+vec_size),0.0,
                           [](const blocked_range<double*>& r, double value)->double {
                               return accumulate(r.begin(),r.end(),value);
                           },plus<double>())/vec_size;
}

double parallelStandardDeviation(double *a, double *meanValue) {
    // reduce by pointing to a place in memory for i in a[] then remove the meanValue and then square
    // emnumerate the above for all a[] then squareroot and divide by array size
    // This is the Lambda syntax,
    return sqrt(parallel_reduce(blocked_range<double*>(a, a+vec_size),0.0,
                           [&](const blocked_range<double*>& r, double sum)->double {
                               for(double* i=r.begin(); i!=r.end(); ++i)
                                   sum += pow(*i-*meanValue,2.0);
                               return sum;
                           }, [](double x, double y)->double {return x+y;}
    )/vec_size);
}

double parallelPearsonCorrelationCoefficient(double *a, double *b) {
    
    double meanA = parallelMean(a);
    double meanB = parallelMean(b);
    
    double standardDeviationA = parallelStandardDeviation(a,&meanA);
    double standardDeviationB = parallelStandardDeviation(b,&meanB);
    // TOP of fraction
    double value = parallel_reduce(blocked_range<size_t>(0, vec_size, 10000), double(0),
                                   [=](blocked_range<size_t> &r, double sum) -> double {
                                       for(size_t i=r.begin();i!=r.end();++i){
                                           sum += ((a[i]-meanA) * (b[i]-meanB));
                                       }
                                       return sum;
                                   },
                                   [=](double a, double b){
                                       return a+b;
                                   });
    
    value *= 1.0/vec_size;
    // BOTTOM of fraction
	value /= (standardDeviationA*standardDeviationB);
    
    return value;
}

double mean(double *a){
	double value = 0.0;
	for(size_t i = 0; i < vec_size; i++) {
		value += a[i];
	}
	return (value/vec_size);
}

double standardDeviation(double *a, double *meanValue){
	double value = 0.0;
	for(size_t i =0; i < vec_size; i++) {
		value += pow((a[i]-*meanValue), 2.0);
	}
	return sqrt(value/vec_size);
}

double pearsonCorrelationCoefficient(double *a, double *b) {
    
    double meanA = mean(a);
    double meanB = mean(b);
    
    double standardDeviationA = standardDeviation(a,&meanA);
    double standardDeviationB = standardDeviation(b,&meanB);
    
    // TOP of fraction
    double value = 0.0;
	for(size_t i =0; i < vec_size; i++) {
		value += ((a[i] - meanA) * (b[i] - meanB));
	}
	value *= 1.0/vec_size;
    
    // BOTTOM of fraction
	value /= (standardDeviationA*standardDeviationB);
    
    return value;
}