import sys
import csv
import operator
import datetime 

max=0
min=9999
f1=csv.reader(open(+str(sys.argv[1])),delimiter=',')

for line in f1:
	numcols=len(line)
	for i in range(0,numcols):
		number=float(line[i])
		if number>max:
			max=number
		if number<min:
			min=number

f1=csv.reader(open(str(sys.argv[1])),delimiter=',')
f2=open(str(sys.argv[2]), 'w')
for line in f1:
	numcols=len(line)
	for i in range(0,numcols):
		number=float(line[i])
		number=int(round((((number-min)/(max-min))*63),0))+1
		f2.write(str(number))
		if(i==(numcols-1)):
			f2.write("\n")
		else:
			f2.write(",")

