#!/bin/bash
# usage: resume.sh REP RESTARTFILE
cd /home/claude/md/rep$1
cp ../in.resume .
echo "resume $2 $(date +%s)" >> resumes.txt
lmp -in in.resume -var rfile $2 -log log.resume.$(date +%s) -screen none
date +%s > t_end
