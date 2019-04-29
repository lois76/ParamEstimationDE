from subprocess import run, Popen, PIPE
import multiprocessing
import time
import csv


def nbContainers():
	result=run("docker ps -aq | wc -l", shell=True, stdout=PIPE).stdout.decode('utf-8')
	return int(result)-1

scilabScriptName = "VoltageClampRIM.sce"
nbCpus = multiprocessing.cpu_count()
maxNbSimPerParam = 20
simId = 0
with open('params.csv') as csv_file:
	csv_reader = csv.reader(csv_file, delimiter=',')
	line_count = 0
	for row in csv_reader:
		# First line contains the header
		if line_count == 0:
			line_count += 1
		else:
			nbSimPerParam = 0
			while nbSimPerParam < maxNbSimPerParam:
				resultPath = "Result_"+row[0]+"_"+str(nbSimPerParam)+"_"+row[1]+"_"+row[2]+"_"+row[3]+"_"+row[4]
				run("mkdir Results/"+resultPath, shell = True)
				command = ["docker run --name scilab"+str(simId)+" --rm  -v $(pwd)/scilab-scripts:/scilab-scripts -v $(pwd)/Results/"+resultPath+":/Result scilab "+scilabScriptName+" "+row[1]+" "+row[2]+" "+row[3]+" "+row[4] + " > outputs/output_docker_"+str(simId)]
				Popen(command, shell = True)
				simId += 1
				nbSimPerParam += 1

				# Loop if the current number of running containers is greater or equal than the number of available cpus
				while nbCpus <= nbContainers():
					time.sleep(1)

