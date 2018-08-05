import os
'''
leetcode 3Sum Closest
2018.8.3
'''
def threesumclosest(nums, target):
    """
    :type nums: List[int]
    :type target: int
    :rtype: int
    """
    nums.sort()
    lenth = len(nums)
    tmp = 0
    if lenth<3:
        return tmp
    com = float("inf")
    for i in range(0,lenth-2):
        l = i+1
        r = lenth-1
        while l<r:
            temp = nums[i]+nums[r]+nums[l]
            var = abs(temp-target)
            if var < com:
                com = var
                tmp = temp
            if temp<target:
                l = l+1
                continue
            elif temp == target:
                break
            else:
                r = r-1
                continue

    return tmp

def help():
    nums = [0,1,2]
    target = 0
    res = threesumclosest(nums,target)
    print res

help()

#if __name__ == "__main__":
#    help()