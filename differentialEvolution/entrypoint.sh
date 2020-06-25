#!/bin/bash

# $2=nbPop, $3=nbIter, $4=F, $5=CR

scilab-adv-cli -f main.sce -nb -args $1 $2 $3 $4 $5 #first arg is the name of script containing the simulation, other args are the parameters 


