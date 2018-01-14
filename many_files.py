from os import walk

import time

start_time = time.time()
c1 = 0
with open('./big_file.txt', 'r') as f:
    for line in f.readlines():
        c1 += 1

print("big file: %d lines, %s seconds" % (c1, time.time() - start_time))

start_time = time.time()
c2 = 0
for (dirpath, dirnames, filenames) in walk('./many_files/'):
    for file in filenames:
        with open('./many_files/' + file) as f:
            for line in f.readlines():
                c2 += 1

print("many files: %d lines, %s seconds" % (c2, time.time() - start_time))
