import csv
import matplotlib
import matplotlib.pyplot as plt

#Parse the scaling result csv file
ux = []  #x-coordinates
uy1 = []  #uniform 1st derivative
uy2 = []  #uniform 2nd derivative
gx = []
gy1 = []  #gaussian 1st derivative
gy2 = []  #gaussian 2nd derivative
modeCount = 0;

with open('scalingTestResult.csv', mode = 'r') as infile:
  csvFile = csv.DictReader(infile, fieldnames = ['mode', 'numSplines', 'paraCoords', '1stDerivDuration', '2ndDerivDuration'])
  for row in csvFile:
    if row['mode'] == 'uniform':
      ux.append(int(row['numSplines']))
      uy1.append(float(row['paraCoords'])*int(row['numSplines'])/float(row['1stDerivDuration']))
      uy2.append(float(row['paraCoords'])*int(row['numSplines'])/float(row['2ndDerivDuration']))
    else:
      gx.append(int(row['numSplines']))
      gy1.append(float(row['paraCoords'])*int(row['numSplines'])/float(row['1stDerivDuration']))
      gy2.append(float(row['paraCoords'])*int(row['numSplines'])/float(row['2ndDerivDuration']))

#Plotting 1st deriv result uniform
fig, ax = plt.subplots()
ax.plot(ux, uy1, 'o', linestyle='-')
ax.set_xlabel("number of splines")
ax.set_ylabel("number of para coords evaluated per second")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(bottom=1)
ax.set_ylim(top=1e15)
ax.set_title("1st Derivative Result with Uniform Spline Length")
offset = 8
for i in range(len(uy1)):
  ax.annotate(
    f"{round(uy1[i], 2)}",
    (ux[i], uy1[i]),
    xytext=(0, offset),
    textcoords='offset points',
    ha='center',
    va='bottom',
    fontsize=10)
  if i % 2 == 0:
    offset = -16
  else:
    offset = 8
plt.show()
fig.savefig('1stDerivUniform.png')

#Plotting 2nd deriv result uniform
fig, ax = plt.subplots()
ax.plot(ux, uy2, 'o', linestyle='-')
ax.set_xlabel("number of splines")
ax.set_ylabel("number of para coords evaluated per second")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(bottom=1)
ax.set_ylim(top=1e15)
ax.set_title("2nd Derivative Result with Uniform Spline Length")
offset = 8
for i in range(len(uy2)):
  ax.annotate(
    f"{round(uy2[i], 2)}",
    (ux[i], uy2[i]),
    xytext=(0, offset),
    textcoords='offset points',
    ha='center',
    va='bottom',
    fontsize=10
  )
  if i % 2 == 0:
    offset = -16
  else:
    offset = 8
plt.show()
fig.savefig('2ndDerivUniform.png')

#Plotting 1st deriv result gaussian
fig, ax = plt.subplots()
ax.plot(gx, gy1, 'o', linestyle='-')
ax.set_xlabel('number of splines')
ax.set_ylabel('number of para coords evaluated per second')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(bottom=1)
ax.set_ylim(top=1e15)
ax.set_title("1st Derivative Result with Variable Spline Length Sampled from a Gaussian Distribution")
offset = 8
for i in range(len(gy1)):
  ax.annotate(
    f"{round(gy1[i], 2)}",
    (gx[i], gy1[i]),
    xytext=(0, offset),
    textcoords='offset points', 
    ha='center',
    va='bottom',
    fontsize=10
  )
  if i % 2 == 0:
    offset = -16
  else:
    offset = 8
plt.show()
fig.savefig('1stDerivGaussian.png')

#Plotting 2nd deriv result gaussian
fig, ax = plt.subplots()
ax.plot(gx, gy2, 'o', linestyle='-')
ax.set_xlabel("number of splines")
ax.set_ylabel("number of para coords evaluated per second")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(bottom=1)
ax.set_ylim(top=1e15)
ax.set_title("2nd Derivative Result with Variable Spline Length Sampled from a Gaussian Distribution")
offset = 8
for i in range(len(gy2)):
  ax.annotate(
    f"{round(gy2[i], 2)}",
    (gx[i], gy2[i]),
    xytext=(0,offset),
    textcoords='offset points',
    ha='center',
    va='bottom',
    fontsize=10)
  if i % 2 == 0:
    offset = -16
  else:
    offset = 8
plt.show()
fig.savefig('2ndDerivGaussian.png')
