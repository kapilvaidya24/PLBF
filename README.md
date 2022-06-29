# PLBF
Code for ICLR 2020 paper [Partitioned Learned Bloom Filter(PLBF)](https://arxiv.org/abs/2006.03176).


plbf_algo.py: Contains the algorithm to find the optimal parameter values.

Input:
- Score samples for keys and non-keys
- Number of partitions in PLBF (5-10 recommended)
- Target False Positive Rate 

Output:
- Score Thresholds 
- False positive rate for each partition

Eg,
Lets suppose we get the following output from the algorithm for K=5.

Thresholds:[0.2, 0.4, 0.6, 0.8, 1.0]
FPR: [0.01, 0.2, 0.5, 1, 1]

This means that inputs(s) with score s(x)<= 0.2 go to the 1st parition, 0.2< s(x) <= 0.4 to the 2nd parition and so on.
The FPR for 1st partition is  0.01, 2nd parition is 0.2 and so on
