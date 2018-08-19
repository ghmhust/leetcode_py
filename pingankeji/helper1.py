def lengthOfLastWord(s):
    if ' ' not in s:
        return len(s)
    A = s.split(' ')
    print A
    while '' in A:
        A.remove('')
    var = len(A)
    if var == 1:
        return len(A[0])
    elif var == 0:
        return 0
    else:
        return len(A[-1])

def main():
    while(1):
        s = raw_input()
        if s == '':
            break
        res = lengthOfLastWord(s)
        print res
    return 0

main()