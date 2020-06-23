from subprocess import run, Popen, PIPE
import multiprocessing
import time
from datetime import datetime
from datetime import timedelta

def nbContainers():
	result=run("docker ps -aq | wc -l", shell=True, stdout=PIPE).stdout.decode('utf-8')
	return int(result)

# Name of the script that should be run
scilabScriptName = "EstimationVoltageRIM11.sce"
nbPop = 190
nbIter = 1000
#mutation = 0.5
#crossover = 0.85

# Starting date
startingDate = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")

# The number of available CPU
nbCpus = multiprocessing.cpu_count()

# The date of the last backup done
lastBackup = datetime.now()

# The index of the first simulation
startingIndex = 1

# The index of the last simulation
endingIndex = 30

# The index of the simulation to start next
nbSimStarted = startingIndex

while nbSimStarted <= endingIndex:

	resultPath = "run_"+str(scilabScriptName)+"_"+str(nbSimStarted)
	run("mkdir Results/"+resultPath, shell = True)

	#command = ["docker run --name scilab"+str(nbSimStarted)+" --rm  -v $(pwd)/scilab-scripts:/scilab-scripts -v $(pwd)/Results/"+resultPath+":/Result scilab "+scilabScriptName+" "+str(nbPop)+" "+str(nbIter)+" "+str(mutation)+" "+str(crossover)+ " > outputs/output_docker_" + str(nbSimStarted) + " 2> outputs/output_docker_error_" + str(nbSimStarted)]
	command = ["docker run --name scilab"+str(nbSimStarted)+" --rm  -v $(pwd)/scilab-scripts:/scilab-scripts -v $(pwd)/Results/"+resultPath+":/Result scilab "+scilabScriptName+" "+str(nbPop)+" "+str(nbIter)+ " > outputs/output_docker_" + str(nbSimStarted) + " 2> outputs/output_docker_error_" + str(nbSimStarted)]
	Popen(command, shell = True)
	nbSimStarted += 1
	time.sleep(30) # wait 30 sec to make sure the container is really started

	# Loop if the current number of running containers is greater or equal than the number of available cpus
	while nbCpus <= nbContainers():
		# Sleep a bit: don't need to be awaken all the time
		time.sleep(350)

while nbContainers() > 0:
	# We started all the containers, but we still have to wait until their end
	time.sleep(350)
