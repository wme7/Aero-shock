def to30(fname):
    inf = open(fname)
    content = inf.read()
    inf.close()
    if '__future__' in content:
        lines = content.split('\n')
        for i, line in enumerate(lines):
            if '__future__' in line:
                lines.pop(i)
                lines.pop(i) # input = raw_input
                if i > 0 and not lines[i-1].strip():
                    lines.pop(i-1)
                outf = open(fname, 'w')
                outf.write('\n'.join(lines))
                outf.close()
                print('Altered', fname)
                return
    else:
        print('No change to', fname)

import sys
if len(sys.argv) > 1 and sys.argv[1].endswith(('.py', '.cgi')):
    to30(sys.argv[1])
else:
    print('Need 2.6 source file parameter!')

            
            
