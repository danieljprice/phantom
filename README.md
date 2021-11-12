Phantom Wind
============

Basic git instructions
----------------------

```
local commit 

- git status
- git diff src/main/eos.F90
- git add filename
- git commit -m "these are all the changes I made made"
- git push

getting a specific revision
- git checkout 536aa21 -b Krome

undo last commit while keeping the changes
- git reset --soft HEAD~1

getting updates:
- git stash
- git pull
- git stash pop

set OMP THREAD number
-  export OMP_NUM_THREADS=1
```

Merge official PHANTOM branch with master
-----------------------------------------
```
- git checkout master
- git pull
- git fetch upstream
- git merge upstream/master
```

Merge local master into wind branch
-----------------------------------
```
- git checkout master
- git pull
- git checkout wind
- git merge master
```
