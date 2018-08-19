import copy
def helper(nums,n,target,temp,res):
    lenth = len(nums)
    if lenth<n:
        return
    l = 0
    r = lenth-1
    if n == 2:
        while l<r:
            if nums[l]+nums[r]<target:
                l = l+1
            elif nums[l]+nums[r]>target:
                r = r-1
            else:
                temp.append(nums[l])
                temp.append(nums[r])
                temp1 = copy.copy(temp)
                res.append(temp1)
                temp.pop()
                temp.pop()
                l = l + 1
        return
    for i in range(0,lenth-1):
        temp.append(nums[i])
        helper(nums[i+1:],n-1,target-nums[i],temp,res)
        temp.pop()

def foursum(nums,target):
    nums.sort()
    lenth = len(nums)
    res = []
    if lenth<4:
        return res
    temp = []
    helper(nums,4,target,temp,res)
    tmp = []
    for i in res:
        tmp.append(str(i))
    tmp = list(set(tmp))
    tmp3 = []
    for i in tmp:
        tmp2 = []
        i = i[1:-1]
        tmp1 = i.split(",")
        for j in tmp1:
            tmp2.append(int(j))
        tmp3.append(tmp2)
    return tmp3

def main():
    a=[0,7,8,8,-10,1,8,-1,7,-10,-7]
    target=7
    res = foursum(a,target)
    print res

main()