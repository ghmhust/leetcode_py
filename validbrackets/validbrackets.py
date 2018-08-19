import copy

def DFS_helper(res,temp,left_size,right_size):
    if left_size == 0 and right_size == 0:
        temp1=copy.copy(temp)
        res.append(temp1)
        return

    if left_size>0:
        DFS_helper(res,temp+'(',left_size-1,right_size)
    if right_size > left_size:
        DFS_helper(res,temp+')', left_size, right_size-1)

def BFS_helper(res,temp,left_size,right_size):
    if left_size == 0 and right_size == 0:
        temp1=copy.copy(temp)
        res.append(temp1)
        return

    if right_size > left_size:
        BFS_helper(res, temp + ')', left_size, right_size - 1)
    if left_size > 0:
        BFS_helper(res, temp + '(', left_size - 1, right_size)

def generateParenthesis(n):
    res = []
    if n == 0:
        return res
    temp = ""
    BFS_helper(res,temp,n,n)
    return res

def main():
    res = generateParenthesis(3)
    print res

main()