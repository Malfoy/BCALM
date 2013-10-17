import sys
with open(sys.argv[1]+".input.dot",'w') as r:
    with open(sys.argv[1]) as f:
        line = f.read()
        line = ''.join([i for i in line if not i.isdigit()])
        for seq in line.split(';'):
            if len(seq) == 0:
                continue
            r.write("%s;\n" % seq)
