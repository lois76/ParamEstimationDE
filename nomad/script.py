from subprocess import run
import csv
import sys

fileName = sys.argv[1]
strArgs = ""
output = ""
with open(fileName) as fileParam:
	reader = csv.reader(fileParam, delimiter=' ')
	for i in next(reader):
		strArgs = strArgs + i + " "
		output = output + "0 "

run("mkdir Result/", shell = True)
command = ["docker run --name scilab_nomad --rm  -v $(pwd)/Result/:/Result scilab_nomad "+strArgs+ " > /dev/null 2>&1"]
run(command, shell = True)
run("docker wait scilab_nomad", shell = True)

with open('Result/result.csv') as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=' ')
	row = next(csv_redear)
	print(row[0] + output)

run("rm -f Result/", shell = True)
