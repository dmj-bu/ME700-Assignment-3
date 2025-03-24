import numpy as np
from setupexample import example_functions as ef


# this is a trivial example of what it may look like to create a tutorial
# you may want to do something like show how to run a function over many examples...
num = 100
x_vec = np.random.random((num, 1))
y_vec = np.random.random((num, 1))
sum_vec = np.zeros((num, 1))
for kk in range(0, num):
    x = x_vec[kk]
    y = y_vec[kk]
    sum = ef.add_x_y(x, y)
    sum_vec[kk] = sum
