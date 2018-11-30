# Expectation Maximization Algorithm

In this code, I have implemented EM algorithm for Bayesian Networks. This is a generalized code which takes data and network(in bif format) and applies EM to learn the CPT.

Here a Bayesian network is provided in bif (Bayesian Interchange Format) which represents a medical diagnosis system. We were also provided with some data records of medical cases (symptoms & diseases), but somehow some of the data was missing, hence EM algorithms comes into the picture. You may need to adapt the data pre-processing code (in diagnoser.cpp) as per the format of your data file.

The main code consists of iterations which terminates when change in probabilities goes below certain threshold (0.0001) in my case. Even in the worst case when CPTs don't converge, code will terminate before 10 minutes, basically a time constraint is used within the loop.

In the starting, all CPT's are initialized randomly using pseudorandom number generator available in c++. Current Time is provided as seed to the rand() so that each time code runs with different initializations and it is really fascinating to see that initialization doesn't matter at all that much (that's the beauty of EM algorithm)
Sometimes it may takes more iteration to converge or score fluctuates a little bit, but it converges to almost same CPT values.

Each iteration consistes of two steps, first one which calculates expectation of the missing data in the given records and the second step consists of maximization step, which basically updates CPT values, by the the prinicple of Maximum Likelihood Estimation, based on the records.

To avoid Nan, zeros and ones in the CPT -> m-smoothing has been used with m = 0.001.

Code has been optimized at various stages, most of the work is done in preprocessing stage (that is before the start of iterations) so that minimal recomputations needed to be done. Currently each iteration takes about 0.3 seconds and preprocessing stage takes about 1s.

For the given data, convergence is attained in about 10 iterations and score obtained is about 19.1 approx.

To run the code, just open the terminal and execute run.sh
To run on your own network and on your own data, just change the files in run.sh and you are done.