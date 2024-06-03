import numpy as np
import random



# binary spliting
def binary_splitting_round(s):
    # s: np.array the infectious status & test status
    num = 0
    flag = sum(s[:,0])>0
    assert flag
    stages = 0
    if len(s[:,0])==1:
        s[0,1] = s[0,0]
        return num,s,stages
    
    B1, B2 = np.array_split(s.copy(), 2,axis=0)
    flag = sum(B1[:,0])>0
    num+=1
    stages += 1
    
    if flag:
        n,stmp,stage = binary_splitting_round(B1)
        s[:len(B1),1] = stmp[:,1]
    else:
        s[:len(B1),1] = 0
        n,stmp,stage = binary_splitting_round(B2)
        s[len(B1):,1] = stmp[:,1]
    num += n
    stages += stage
    return num,s,stages 

def binary_splitting(s):
    # modified bs
    # s: 1-d array the infectious status
    st = np.zeros((len(s),2))
    st[:,0] = s
    st[:,1] = np.nan
    nums = 0
    count = sum(np.isnan(st[:,1]))
    stages = 0
    # the undetermined people
    while count!=0:
        mask = np.isnan(st[:,1])
        flag = sum(st[mask,0]>0)>0
        nums += 1
        stages+=1
        if not flag:
            st[mask,1] = 0
        else:
            n,stmp,stage = binary_splitting_round(st[mask,:])
            st[mask,1] = stmp[:,1]
            nums += n
            stages += stage
        count = sum(np.isnan(st[:,1]))
        
    assert sum(st[:,0]!=st[:,1])==0
    return nums,stages, st[:,1]

# diag
def diagalg_iter(s):
    # s(np.array): binary string of infection status
    k = int(np.log2(len(s)))
    l = int(2**(k-1))
    lp = 0
    p = np.zeros(k+1)
    group = dict()
    num = np.ones(k+1,dtype=np.int32)
    for i in range(k):
        p[i] = sum(s[lp:lp+l])>0
        group[i] = s[lp:lp+l]
        num[i] = l
        lp+=l
        l = l//2

    p[-1] = s[-1]
    group[k] = np.array([s[-1]])
    # p(array): pattern
    # group(dict): indicate the group information
    # num(array): the group size
    return p.astype(np.int32), group,num


def diag_splitting(s):
    # s(np.array): binary string of infection status
    num_tests = 0
    stages = 0
    pattern, group, nums = diagalg_iter(s)
    stages +=1
    num_tests += len(pattern)
    indices = np.where(pattern == 1)[0]
    flag = 0
    for i in indices:
        if nums[i]>1:
            num_test,stage = diag_splitting(group[i])
            num_tests += num_test
            if not flag:
                stages+=stage
                flag = 1
    return num_tests,stages




def Qtesting1(s):
    '''
    s(np.array): binary string of infection status
    '''
    import math
    def helper(s, stage=0):
        total = len(s)
        if total == 1:
            return 1, stage  # Base case: no more tests needed, return current stage
        if sum(s) == total or 0:
            return 1, stage
        if sum(s) == 1 or sum(s) == total - 1:
            return math.log2(total) + 1, stage

        infected = sum(s)
        
        m = max(infected, total - infected)
        if m/total < 3/256:
            m = infected

        x = m / total
        num_tests = 1

        # Find the smallest i such that x <= (2^i - 1) / y
        for i in range(1, int(math.log2(total)) + 1):
            y = 2 ** i
            if x >= (y - 1) / y:
                continue
            else:
                # Split the array into y subgroups
                subgroup_size = total // y
                test_count = 0
                max_stage = 0
                # Recursively calculate the number of tests for each subgroup
                for j in range(0, total, subgroup_size):
                    subgroup_tests, subgroup_stage = helper(s[j:j + subgroup_size], stage + 1)
                    test_count += subgroup_tests
                    max_stage += subgroup_stage
                return num_tests + test_count, max_stage

        return num_tests, stage  # Return the total number of tests including this level and current stage

    total_tests, max_stages = helper(s, 0)
    return total_tests, max_stages
    


def Qtesting2(s):
    '''
    s(np.array): binary string of infection status
    '''
    num_tests = 0
    stages = 0
    import math
    
    def helper2(s, stage=0):
        total = len(s)
        if total == 1:
            return 1, stage  # Base case: no more tests needed, return current stage
        
        if sum(s) == 0:
            return 1, stage
        if sum(s) == 1:
            return math.log2(total) + 1, stage
        
        if sum(s) >= 8:
            val = 5
        elif sum(s) >= 4:
            val = 4
        elif sum(s) >= 2:
            val = 3

        if total > 16:
            if val == 5:
                m = total - 8
            if val == 4:
                m = total - 4
            if val == 3:
                m = total - 3
            x = m / total
            num_tests = 1

            # Find the smallest i such that x <= (2^i - 1) / y
            for i in range(1, int(math.log2(total)) + 1):
                y = 2 ** i
                if x >= (y - 1) / y:
                    continue
                else:
                    # Split the array into y subgroups
                    subgroup_size = total // y
                    test_count = 0
                    max_stage = 0
                    # Recursively calculate the number of tests for each subgroup
                    for j in range(0, total, subgroup_size):
                        subgroup_tests, subgroup_stage = helper2(s[j:j + subgroup_size], stage + 1)
                        test_count += subgroup_tests
                        max_stage += subgroup_stage
                    return num_tests + test_count, max_stage



        infected = sum(s)
        m = max(infected, total - infected)
        x = m / total
        num_tests = 1

        # Find the smallest i such that x <= (2^i - 1) / y
        for i in range(1, int(math.log2(total)) + 1):
            y = 2 ** i
            if x >= (y - 1) / y:
                continue
            else:
                # Split the array into y subgroups
                subgroup_size = total // y
                test_count = 0
                max_stage = 0
                # Recursively calculate the number of tests for each subgroup
                for j in range(0, total, subgroup_size):
                    subgroup_tests, subgroup_stage = helper2(s[j:j + subgroup_size], stage + 1)
                    test_count += subgroup_tests
                    max_stage += subgroup_stage
                return num_tests + test_count, max_stage

        return num_tests, stage  

    
    total_tests, max_stages = helper2(s, 0)
    return total_tests, max_stages



def Qtesting1_comm_aware(s,communities):
    '''
    s(np.array): binary string of infection status
    communities(list): the community sizes
    '''
    num_tests = 0
    stages = 0
    ###################################################
    
    community_indices = list(range(len(communities)))

    def search_communities(indices):
        nonlocal num_tests, stages
        if not indices:
            return []
        elif len(indices) == 1:
            # Directly test this single community if only one left
            index = indices[0]
            start = sum(communities[:index])
            end = start + communities[index]
            if np.any(s[start:end]):
                num_tests += 1
                return [(start, end)]
            return []

        mid = len(indices) // 2
        left_indices = indices[:mid]
        right_indices = indices[mid:]

        # Test the pooled communities on left and right
        left_start = sum(communities[:left_indices[0]])
        left_end = left_start + sum(communities[left_indices[0]:left_indices[-1] + 1])
        right_start = sum(communities[:right_indices[0]])
        right_end = right_start + sum(communities[right_indices[0]:right_indices[-1] + 1])

        infected_communities = []
        if np.any(s[left_start:left_end]):
            num_tests += 1
            infected_communities += search_communities(left_indices)
        if np.any(s[right_start:right_end]):
            num_tests += 1
            infected_communities += search_communities(right_indices)

        return infected_communities

    
    # Start the recursive binary search
    infected_communities = search_communities(community_indices)
    
    # Apply Qtesting1 to each detected infected community
    for start, end in infected_communities:
        community_s = s[start:end]
        community_tests, community_stages = Qtesting1(community_s)
        num_tests += community_tests
        stages = max(stages, community_stages)
    
    ###################################################



    return num_tests,stages

def Qtesting2_comm_aware(s,communities):
    '''
    s(np.array): binary string of infection status
    communities(list): the community information
    '''
    num_tests = 0
    stages = 0
    ###################################################
    community_indices = list(range(len(communities)))

    def search_communities(indices):
        nonlocal num_tests, stages
        if not indices:
            return []
        elif len(indices) == 1:
            # Directly test this single community if only one left
            index = indices[0]
            start = sum(communities[:index])
            end = start + communities[index]
            if np.any(s[start:end]):
                num_tests += 1
                return [(start, end)]
            return []

        mid = len(indices) // 2
        left_indices = indices[:mid]
        right_indices = indices[mid:]

        # Test the pooled communities on left and right
        left_start = sum(communities[:left_indices[0]])
        left_end = left_start + sum(communities[left_indices[0]:left_indices[-1] + 1])
        right_start = sum(communities[:right_indices[0]])
        right_end = right_start + sum(communities[right_indices[0]:right_indices[-1] + 1])

        infected_communities = []
        if np.any(s[left_start:left_end]):
            num_tests += 1
            infected_communities += search_communities(left_indices)
        if np.any(s[right_start:right_end]):
            num_tests += 1
            infected_communities += search_communities(right_indices)

        return infected_communities

    
    # Start the recursive binary search
    infected_communities = search_communities(community_indices)
    
    # Apply Qtesting2 to each detected infected community
    for start, end in infected_communities:
        community_s = s[start:end]
        community_tests, community_stages = Qtesting2(community_s)
        num_tests += community_tests
        stages = max(stages, community_stages)
    ###################################################



    return num_tests,stages