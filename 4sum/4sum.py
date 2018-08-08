import copy
def helper(nums,n,target,temp,res,c,dic):
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
                #之前是append的temp 出错
                c1 = copy.copy(c)
                c1 = c1 + str(nums[l])+"#"+str(nums[r])+"#"
                if c1 not in dic:
                    dic[c1] = temp1
                    res.append(temp1)
                temp.pop()
                temp.pop()
                l = l + 1
        return
    for i in range(0,lenth-1):
        temp.append(nums[i])
        c1 = copy.copy(c)
        c = c+str(nums[i])+"#"
        helper(nums[i+1:],n-1,target-nums[i],temp,res,c,dic)
        temp.pop()
        c = copy.copy(c1)

def foursum(nums,target):
    nums.sort()
    lenth = len(nums)
    res = []
    if lenth<4:
        return res
    temp = []
    dic = {}
    c = ""
    helper(nums,4,target,temp,res,c,dic)
    print dic
    return res

def main():
    a=[0,7,8,8,-10,1,8,-1,7,-10,-7]
    target=7
    res = foursum(a,target)
    print res

main()