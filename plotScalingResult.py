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
      uy1.append(float(row['paraCoords'])/float(row['1stDerivDuration']))
      uy2.append(float(row['paraCoords'])/float(row['2ndDerivDuration']))
    else:
      gx.append(int(row['numSplines']))
      gy1.append(float(row['paraCoords'])/float(row['1stDerivDuration']))
      gy2.append(float(row['paraCoords'])/float(row['2ndDerivDuration']))

#Plotting 1st deriv result uniform
plt.plot(ux, uy1, 'o', linestyle='-')
for i in range(len(uy1)):
  label= f"({ux[i]}, {round(uy1[i], 3)})"
  plt.text(ux[i], uy1[i]+3, label, ha='center', va='bottom', fontsize=5)
plt.xlabel("number of splines")
plt.ylabel("number of para coords evaluated per second")
plt.title("1st Derivative Result with Uniform Spline Length")
plt.yscale('log')
plt.ylim(bottom=1)
plt.ylim(top=1000000)
plt.show()
plt.savefig('1stDerivUniform.png')

#Plotting 2nd deriv result uniform
plt.plot(ux, uy2, 'o', linestyle='-')
for i in range(len(uy2)):
  label=f"({ux[i]}, {round(uy2[i], 3)})"
  plt.text(ux[i], uy2[i]+3, label, ha='center', va='bottom', fontsize=5)
plt.xlabel("number of splines")
plt.ylabel("number of para coords evaluated per second")
plt.title("2nd Derivative Result with Uniform Spline Length")
plt.yscale('log')
plt.ylim(bottom=1)
plt.ylim(top=1000000)
plt.show()
plt.savefig('2ndDerivUniform.png')

#Plotting 1st deriv result gaussian
plt.plot(gx, gy1, 'o', linestyle='-')
for i in range(len(gy1)):
  label=f"({gx[i]}, {round(gy1[i], 3)})"
  plt.text(gx[i], gy1[i]+3, label, ha='center', va='bottom', fontsize=5)
plt.xlabel("number of splines")
plt.ylabel("number of para coords evaluated per second")
plt.title("1st Derivative Result with Variable Spline Length Sampled from a Gaussian Distribution")
plt.yscale('log')
plt.ylim(bottom=1)
plt.ylim(top=1000000)
plt.show()
plt.savefig('1stDerivGaussian.png')

#Plotting 2nd deriv result gaussian
plt.plot(gx, gy2, 'o', linestyle='-')
for i in range(len(gy2)):
  label=f"({gx[i]}, {round(gy2[i], 4)})"
  plt.text(gx[i], gy2[i]+3, label, ha='center', va='bottom', fontsize=5)
plt.xlabel("number of splines")
plt.ylabel("number of para coords evaluated per second")
plt.title("2nd Derivative Result with Variable Spline Length Sampled from a Gaussian Distribution")
plt.yscale('log')
plt.ylim(bottom=1)
plt.ylim(top=1000000)
plt.show()
plt.savefig('2ndDerivGaussian.png')
