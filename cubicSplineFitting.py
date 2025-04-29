import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

import csv
x = []
y = []
with open('inputPts.csv', mode ='r') as file:    
  csvFile = csv.DictReader(file, fieldnames = ['x', 'y', 'z'])
  for lines in csvFile:
    x.append(float(lines['x']))
    y.append(float(lines['y']))
    print(lines)


csc = CubicSpline(x, y, bc_type='clamped')

print(f"x: {x}")
print(f"y: {y}")
print(f"yInterp: {csc(x)}")

fig, ax = plt.subplots(figsize=(8, 10))
ax.plot(x, y, 'o', label='data')
ax.plot(x, csc(x), label="clamped")
ax.legend(loc='lower left', ncol=1)
plt.show()
plt.savefig('cubicSpline.png')
