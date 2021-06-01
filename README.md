Phantom AGB
===========


**+++++ Git instructions +++++**

**1) Commit locally:**
--------------------

- git status
- git diff
- git add .
- git commit -m "<comment>"

**2) Push to remote repository:**
-----------------------------

- git push phantom AGBlab

**3) Sync branch with original master:**
-------------------------------------

- git checkout master
- git fetch upstream
- git merge upstream/master

  **4) Sync branch with master:**
  ----------------------------
  
- git checkout master
- git pull
- git checkout AGBlab
- git merge master
