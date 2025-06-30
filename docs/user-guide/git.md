# Working with git and MOCCa

Disclaimer: this guide is solely intended to get people without git experience started.

Make sure that 

1. your local pc has git installed
2. have a github account. 

In my experience, setting up ssh as your method of access to github is a good idea; see [here](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account) for more information.

## Getting a brand new copy of MOCCa

Execute the following command in a directory of your choice

    git clone https://github.com/IAA-nuclear/tantalus_full.git

If your computer succesfully communicates with github, you will now have a ```tantalus_full``` 
directory that contains the code.

## Setting your username and email address

To sign off on your work, git needs to know who you are:

    git config --global user.name "Jan Jaak"
    git config --global user.email "jan.jaak@ulb.be"

## Receiving updates from your current branch

When you want to receive an update from the repository, the recommended 
workflow is the following:

     git stash
     git pull
     git stash pop

The first command - ``` stash ``` - temporarily removes your (uncommitted) local changes. 
The second will ask github for updates. ```stash pop``` will reapply the local changes
you stored in the first step. 
          
## Checking out a particular MOCCa version 

If you want to work with a specific version of the code, you can use the following commands

    git fetch 
    git checkout vA.B.C

The first command updates your local pc with all relevant tags; the second ensures 
that your local copy of the code will be ``` vA.B.C```.

You can find the release notes and the logic behind the versioning tags on [this page](../release.md).
