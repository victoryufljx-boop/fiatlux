#!/bin/bash
# usage: final.sh [fast]  -> final assemble + push
set -e
cd /home/claude/md
O=/home/claude/fiatlux/md_sites
if [ "$1" = fast ]; then
  nice python3 analyze.py anafast.json rep5 rep6 > anafast.txt
  FAST="rep5 rep6" FASTANA=anafast.json python3 assemble.py ana1234.json $O/Er-Ge10-SiO2.json rep1 rep2 rep3 rep4
  python3 toxyz.py rep5 rep6
else
  python3 assemble.py ana1234.json $O/Er-Ge10-SiO2.json rep1 rep2 rep3 rep4
fi
python3 make_md.py $O/Er-Ge10-SiO2.json $O/Er-Ge10-SiO2.md
cp *.py in.* *.sh $O/scripts/
cd /home/claude/fiatlux
git add md_sites
git commit -q -m "md2: Er-Ge10-SiO2 final sites and write-up (Pedone 2006, 4 replicas at 10 K/ps${1:+ + 2 fast-quench replicas})

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01GmgwxcZpDZuJmGJVMUD4qP"
git push -q
git log --oneline -1
