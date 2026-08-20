import csv
import matplotlib.pyplot as plt

#Parse the scaling result csv file
x = []
y1 = []
y2 = []
modeCount = 0;

with open('scalingTestResult.csv', mode = 'r') as infile:
  csvFile = csv.DictReader(infile, fieldnames = ['mode', 'numSplines', 'paraCoords', '1stDerivDuration', '2ndDerivDuration'])
  for row in csvFile:
    x.append(int(row['numSplines']))
    y1.append(float(row['paraCoords'])/float(row['1stDerivDuration']))
    y2.append(float(row['paraCoords'])/float(row['2ndDerivDuration']))
    if (row['mode'] == 'uniform'):
      modeCount += 1

print("x: ", x)
print("y1: ", y1)
print("y2: ", y2)
print("modeCount: ", modeCount)

#Plotting 1st deriv result uniform
plt.plot(x[:modeCount], y1[:modeCount], 'o', linestyle='-')
for i in range(modeCount):
  label= f"{y1[i]}"
  plt.text(x[i], y1[i]+1, label, ha='center', va='bottom', fontsize=10)
plt.xlabel("number of splines")
plt.ylabel("number of para coords evaluated per second")
plt.title("1st Derivative Result with Uniform Spline Length")
plt.show()
plt.savefig('1stDerivUniform.png')

#Plotting 2nd deriv result uniform
plt.plot(x[:modeCount], y2[:modeCount], 'o', linestyle='-')
for i in range(modeCount):
  label=f"{y2[i]}"
  plt.text(x[i], y1[i]+1, label, ha='center', va='bottom', fontsize=10)
plt.xlabel("number of splines")
plt.ylabel("number of para coords evaluated per second")
plt.title("2nd Derivative Result with Uniform Spline Length")
plt.show()
plt.savefig('2ndDerivUniform.png')

#Plotting 1st deriv result gaussian
plt.plot(x[:modeCount], y1[modeCount: len(y1)], 'o', linestyle='-')
for i in range(modeCount, len(y1)):
  label=f"{y1[i]}"
  plt.text(x[i], y1[i]+1, label, ha='center', va='bottom', fontsize=10)
plt.xlabel("number of splines")
plt.ylabel("number of para coords evaluated per second")
plt.title("1st Derivative Result with Variable Spline Length Sampled from a Gaussian Distribution")
plt.show()
plt.savefig('1stDerivGaussian.png')

#Plotting 2nd deriv
plt.plot(x[:modeCount], y2[modeCount: len(y2)], 'o', linestyle='-')
for i in range(modeCount, len(y2)):
  label=f"{y2[i]}"
  plt.text(x[i], y2[i]+1, label, ha='center', va='bottom', fontsize=10)
plt.xlabel("number of splines")
plt.ylabel("number of para coords evaluated per second")
plt.title("2nd Derivative Result with Variable Spline Length Sampled from a Gaussian Distribution")
plt.show()
plt.savefig('2ndDerivGaussian.png')


