import subprocess
import os
def main():
    #p = subprocess.Popen(temp)
    out_file = open('/Users/guohaomin/a','w')
    output = os.popen("python /Users/help.py")
    out_file.close()
    #return_code = p.returncode
    #print "aaaa %s" % p.returncode
    return 0

main()