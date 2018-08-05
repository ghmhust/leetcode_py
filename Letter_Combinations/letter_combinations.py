'''
leetcode Letter Combinations of a Phone Number
2018.8.3
'''
def help(call,num,digits,res,temp):
    if num == len(digits):
        res.append(temp)
        return
    index = int(digits[num])-2
    for i in range(0,len(call[index])):
        var = temp
        temp = temp+call[index][i]
        num = num + 1
        help(call,num,digits,res,temp)
        temp = var
        num = num - 1

def lettercombinations(digits):
    """
    :type digits: str
    :rtype: List[str]
    """
    call = ['abc','def','ghi','jkl','mno','pqrs','tuv','wxyz']
    lenth = len(digits)
    res = []
    num = 0
    temp = ''
    help(call,num,digits,res,temp)
    print res

if __name__ == "__main__":
    test = '233'
    lettercombinations(test)