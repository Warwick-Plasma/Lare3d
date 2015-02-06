WORKING WITH THE GIT REPOSITORY
===============================

The easiest method of downloading the code is to grab a copy of the latest
release package, located in the "Releases" section of the CCPForge site.

For more advanced users, the code is also hosted on a git repository. Details
can be found in the "Git" section on the CCPForge site. There is quite a
steep learning curve for using git, so using this repository is only
recommended for more advanced users who are comfortable that they can deal with
a "git conflict".

One other added complication, is that the LareXd repositories also use git
submodules for tracking the SDF file format. This adds an extra source of
possible issues. However, once a workflow is established it can all be quite
straightforward to work with.

To perform an initial checkout of the code using git, you should issue the
following command:

  git clone --recursive ssh://<name>@ccpforge.cse.rl.ac.uk/gitroot/lare3d

The "--recursive" flag ensures that not only the "lare3d"
repository is checked out, but also the "SDF" submodules.

It is recommended that after checking out a copy of the git repository, users
immediately create a new working branch and leave the default "master" branch
untouched. A new branch can be created and switched to with the command
"git checkout -b work".

When you wish to update to the latest version, do the following sequence of
actions. First, commit or stash any changes you have made in your "work"
branch. Next, switch to the "master" branch with
"git checkout master". Now pull the changes with "git pull",
followed by "git submodule update --recursive".
At this stage your "master" branch should be fully up to date.

Merging the new version in with your "work" branch is prone to error, so it
is recommended that you create a temporary copy of this branch just in case
everything goes wrong. The command "git branch workold work" will
create a branch named "workold" which is just a copy of "work". This branch
can be deleted once the merge is completed successfully. If everything goes
wrong in the "work" branch, you can reset it back to the original using the
command "git reset --hard workold".

In order to update your work branch, switch back to it with
"git checkout work" and merge in the changes with "git merge master".
After issuing this last command, there is a fair chance that you will encounter
conflicts. You must now resolve those conflicts and commit the changes.
After successfully merging in the changes, you can now delete the temporary
copy of your work branch with "git branch -D workold".

