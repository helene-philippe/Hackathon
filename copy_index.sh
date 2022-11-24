#!/bin/bash

for file in {"Genome","Log.out","SA","SAindex","chrLength.txt","chrName.txt","chrNameLength.txt","chrStart.txt","genomeParameters.txt"}
	do 
	scp  ubuntu@134.158.249.51:hackathon/ref/$file /Users/helenephilippe/Desktop/3A/AMI2B/hackathon/Hackathon/ref/$file
	done  