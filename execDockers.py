from subprocess import Popen, PIPE
import multiprocessing
import time

scilabScriptName="EstimationSteadyCurrent11.sce"
nbCpus=multiprocessing.cpu_count()
nbContainers=0
nbSim=0
while nbSim<100:
	pipe=Popen("docker container ls | wc -l", shell=True, stdout=PIPE)
	nbContainer=int(pipe.communicate()[0])-1
	if nbContainers<nbCpus :
		Popen("mkdir Result"+str(nbSim), shell=True)
		command=["docker run --name scilab"+str(nbSim)+" --rm -v $(pwd)/Result"+str(nbSim)+":/Result scilab "+scilabScriptName+" 10 10 0.5 0.9"]
		Popen(command, shell=True)
		nbSim+=1
	else:
		time.sleep(100)
