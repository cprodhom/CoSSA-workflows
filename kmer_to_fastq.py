#! /usr/bin/python

import fileinput

nr=0
for line in fileinput.input():
    pass

    if line.startswith(('#','+')):
        continue
    else:
        stripline = line.rstrip()
        splitline = stripline.split()
        head = nr
        read = splitline[0]
        qual = "I" * len(splitline[0])
        print(">{0}\n{1}\n+\n{2}".format(head, read, qual))
        nr += 1

