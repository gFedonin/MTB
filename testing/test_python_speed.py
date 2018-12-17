from datetime import datetime

a = range(10000)

start_time = datetime.now()
std = sum(i*i for i in a)
mean = sum(a)
print('using list comprehension to sum: {}'.format(datetime.now() - start_time))


start_time = datetime.now()
std = 0
mean = 0
for x in a:
    mean += x
    std += x*x
print('using cycle to sum: {}'.format(datetime.now() - start_time))