
#!/usr/bin/env python2.7

import optparse

def main() :
    parser = optparse.OptionParser(version='%prog version 1.0')
    parser.add_option('-f', '--infile', help='input pdb file)', action='store')
    options, arguments = parser.parse_args()

#***********************************************************************************************

    infile =open(str(options.infile),"r")
    ow=[]
    hw1=[]
    hw2=[]
    p=[]; w=[]
    for line in infile:
        if line.strip() :
            fields=line.split()
            if (len(fields)>1) and fields[0]=="ATOM":
                if (fields[2]=="OW"):
                    ow.append(fields[1])
                if (fields[2]=="HW1"):
                    hw1.append(fields[1])
                if (fields[2]=="HW2"):
                    hw2.append(fields[1])
                if (fields[2]=="W"):
                    w.append(fields[1])
    infile.close()
    print('w:',(int(len(w))))
    print ('ow:', (int(len(ow))))
    print ('hw1:', (int(len(hw1))))
    print ('hw2:', (int(len(hw2))))
    print ('lipidos:', (int(len(p))))

if __name__=="__main__" :
    main()
