from subprocess import run
import csv
import sys

strArgs = ""
output = ""
for arg in sys.argv:
	strArgs = strArgs + arg + " "
	output = output + " 0"

run("mkdir Result/", shell = True)
command = ["docker run --name scilab_nomad --rm  -v $(pwd)/Result/:/Result scilab_nomad "+strArgs+ " > /dev/null 2>&1"]
run(command, shell = True)

with open('Result/result.csv') as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=' ')
	print next(csv_reader)[0] + output

