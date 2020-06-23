#!/bin/bash

# $2=taillePop, $3=nbiter, $4=F, $5=CR

scilab-adv-cli -f main.sce -nb -args $1 $2 $3 #first arg is the name of script containing the simulation, other args are the parameters 


