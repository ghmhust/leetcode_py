def main():
    f = ["123213333.request","dawdwcv333.request","dadawdwad","adawxxx33"]
    res = []
    for item in f:
        if len(item)>len(".request"):
            if item[-len(".request"):] == ".request":
                res.append(item)
    print res
    return 0

main()