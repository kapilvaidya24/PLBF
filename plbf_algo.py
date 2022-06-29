

import math
import numpy as np
import operator

#Given the key, non-key scores and num partitions, gives threshold values as output
def get_thresholds(key_scores,non_key_scores,num_partitions):

    #Chunking the score space. The naive DP algorithm is quadratic in size of the set(too expensive)
    complete_list=key_scores+non_key_scores
    complete_list.sort()

    #using 100 chunks or lesser if needed
    # Larger chunks miifght lead to an error be careful
    num_chunks=min(100,int(len(complete_list)/5))
    chunk_thresholds=[0.0]*num_chunks

    
    for i in range(0,num_chunks):
        index=len(complete_list)*1.00*(i+1)/num_chunks
        index=int(min(index,len(complete_list)-1))
        chunk_thresholds[i]=complete_list[index]

    chunk_thresholds=list(set(chunk_thresholds))
    chunk_thresholds.sort()
    num_chunks=len(chunk_thresholds)

    chunk_dict={}
    for i in range(0,num_chunks):
        chunk_dict[i]=chunk_thresholds[i]


    #Computing the key aand non-key density for each chunk
    G_arr=[0.0]*num_chunks
    H_arr=[0.0]*num_chunks
    G_cummulative=[0.0]*num_chunks
    H_cummulative=[0.0]*num_chunks

    for i in range(0,len(key_scores)):
        for j in range(0,len(chunk_thresholds)):
            if key_scores[i] <= chunk_thresholds[j]:
                G_arr[j]+=(1.00/len(key_scores))
                break

    for i in range(0,len(non_key_scores)):
        for j in range(0,len(chunk_thresholds)):
            if non_key_scores[i] <= chunk_thresholds[j]:
                H_arr[j]+=(1.00/len(non_key_scores))
                break       

    for i in range(0,num_chunks):
        if i==0:
            G_cummulative[i]=G_arr[i]            
            H_cummulative[i]=H_arr[i]            
        else:
            G_cummulative[i]=G_arr[i]+G_cummulative[i-1]            
            H_cummulative[i]=H_arr[i]+H_cummulative[i-1]
      
    #DP algorithm starts here. 
    # DP[i][j]: What is the maximum KL divergence while considering i partitions of first j chunks.
    # jump_matrix[i][j]: What was the endpoint of the last partition for maximum KL divergence of DP[i][j]
    DP_matrix=[]
    jump_matrix=[]
    for i in range(0,num_partitions):
        temp=[0.0]*num_chunks
        temp1=[0.0]*num_chunks
        DP_matrix.append(temp)
        jump_matrix.append(temp1)


    for i in range(0,num_chunks):
        DP_matrix[0][i]=G_cummulative[i]*math.log(G_cummulative[i]*1.00/H_cummulative[i])
        jump_matrix[0][i]=0-1


    for i in range(1,num_partitions):
        for j in range(0,num_chunks):
            if j<i:
                DP_matrix[i][j]=0-1
                jump_matrix[i][j]=0-1
                continue

            poss_list={}       
            for k in range(0,j):
                if k <(i-1):
                    continue
                delta_G=G_cummulative[j]-G_cummulative[k]
                delta_H=H_cummulative[j]-H_cummulative[k]    
                ans=DP_matrix[i-1][k]+(delta_G)*math.log((delta_G)*1.00/(delta_H))
                poss_list[k]=ans

            #Find the partition that resulted in maximum KL Divergence
            jump_loc=max(poss_list, key=poss_list.get) 
            
            DP_matrix[i][j]=poss_list[jump_loc]
            jump_matrix[i][j]=jump_loc

    # Using the jump matrix to find which chunks lead to maximum KL divergence. This is used to generate the threshold arrays
    threshold_params=[1.00]
    curr_loc=num_chunks-1
    for i in range(0,num_partitions-1):
        curr_loc=jump_matrix[num_partitions-1-i][curr_loc]
        threshold_params.append(chunk_dict[curr_loc])

    threshold_params.sort()

    print("Max KL Divergence is: ",DP_matrix[num_partitions-1][num_chunks-1])
    print("Threshold params:",threshold_params)


    return threshold_params


#Gets the optimal FPR of regions given the threshold, scores and target false positive rate.
def get_optimal_fpr(key_scores,non_key_scores,threshold_params,target_fpr):


    # Simply calculating density of keys/noon-keys in each region
    threshold_params.sort()
    num_partitions=len(threshold_params)

    G_arr=[0.0]*num_partitions
    H_arr=[0.0]*num_partitions
    fpr_params=[1.0]*num_partitions

    for i in range(0,len(key_scores)):
        for j in range(0,len(threshold_params)):
            if key_scores[i] <= threshold_params[j]:
                G_arr[j]+=(1.00/len(key_scores))
                break

    for i in range(0,len(non_key_scores)):
        for j in range(0,len(threshold_params)):
            if non_key_scores[i] <= threshold_params[j]:
                H_arr[j]+=(1.00/len(non_key_scores))
                break           


    #Using iterative algorithm to find optimal fpr values. Algorithm 1 in the paper 
    G_sum=0
    H_sum=0

    for i in range(0,num_partitions):
        fpr_params[i]= (target_fpr*G_arr[i]*1.00)/H_arr[i]

    is_fpr_larger_than_one=True
    while(is_fpr_larger_than_one):

        G_sum=0
        H_sum=0
    
        for i in range(0,num_partitions):
            if(fpr_params[i]>=1):
                fpr_params[i]=1
                G_sum+=G_arr[i]
                H_sum+=H_arr[i]

        is_fpr_larger_than_one=False
        calc_fpr=0.0
        for i in range(0,num_partitions):
            calc_fpr+=H_arr[i]*fpr_params[i]
            if fpr_params[i]>=1:
                continue

            new_fpr=((target_fpr-H_sum)*G_arr[i])*1.00/((1-G_sum)*H_arr[i])

            if new_fpr>=1 and fpr_params[i]<1:
                is_fpr_larger_than_one=True

            fpr_params[i]=new_fpr


    
    calc_fpr=0.0

    for i in range(0,num_partitions):
        calc_fpr+=H_arr[i]*fpr_params[i]
    
    #Checking if target fpr matches calculated fpr
    print("Target FPR:",target_fpr," PLPBF estimated FPR:",calc_fpr)

    return fpr_params

#finds optimal threshold and fpr values for PLBF configuration
# NOTE: MAKE SURE ALL SCORES ARE BETWEEN 0 AND 1
# key_scores: It is a list of learned model outputs for keys in the set.  Eg,[1.0,0.9,0.8...]
# non_key_scores: It is a list of learned model outputs for elements not in the set. Eg,[1.0,0.9,0.8...]
# num_partitions: It is the number of pparitions for PLBF structure
# target_fpr: It is the target fpr that is desired by the user. 
def get_parameter_vals(key_scores,non_key_scores,target_fpr,num_partitions):

    threshold_params=get_thresholds(key_scores,non_key_scores,num_partitions)
    fpr_params=get_optimal_fpr(key_scores,non_key_scores,threshold_params,target_fpr)

    return threshold_params,fpr_params




#Generated some sample key scores with mean=0.7 and standard deviation=0.2.
# Key scores would be concentrated towards 1
mu, sigma = 0.7, 0.2 
key_scores = np.random.normal(mu, sigma, 100000)
key_scores=[item for item in key_scores if (item >= 0 and item<=1.0)]

#Generated some sample non-key scores with mean=0.3 and standard deviation=0.2
# Non Key scores would be concentrated towards 0
mu, sigma = 0.3, 0.2 
non_key_scores = np.random.normal(mu, sigma, 100000)
non_key_scores=[item for item in non_key_scores if (item >= 0 and item<=1.0)]

# 5 partitions in PLBF
K=5

# Fixing Target FPR
target_fpr=0.1

threshold_params, fpr_params = get_parameter_vals(key_scores,non_key_scores,target_fpr,K)

print("PLBF Recommended Thresholds: ",threshold_params)
print("PLBF Recommended FPR's: ",fpr_params)
