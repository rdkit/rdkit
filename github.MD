Doc for pull request with squash option during merge

% get the last master official RDKit (code before any modification) 
step 1 
in github web site: select an old version of RDkit and create a branch from that commit
in the terminal local repo: git pull the new branch 
in mac git desktop app: click on the button update the from the rdkit/master code than on the sync button

% get the last branch sync with rdkit/master
step 2 in the terminal local repo: git checkout "the new branch"

% merge local (origin) master (mine!) all in one to the "new branch" to have only one commit with all my modifications!
step 3 in the temincal local repo: git merge --squash master

% commit the last merge
step 4 in the temincal local repo: git commit -m "3D Descriptors Dragons"

% push the last commit into the github web:
step 5 in the terminal local repo: git push

this is it!
