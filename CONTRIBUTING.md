# Contributing Guidelines
There are many ways to contribute to MBDyn:

 - using it and reporting bugs, fixing/adding tutorials, documentation, and so on;
 - submitting patches that fix bugs or implement new features, and so on;
 - establishing a grant with DAER/Polimi to implement the features or develop 
   the models you need.

# Access to the Gitlab repository
The MBDyn main repository was been recently moved to this Gitlab repository, 
hosted on the 
[POLIMI](https://www.polimi.it/) servers. At the moment, the access to some
Gitlab functionalities for external users is limited: we're working with 
the ITC department of [POLIMI](https://www.polimi.it/) to improve it. 
However, in the meantime, it is **necessary,
even for registered POLIMI or Google users, to be registered as Project Members**
for them to be able to open issues, fork the project and therefore also to submit
merge requests. Please contact us through the 
[MBDyn users mailing list](https://www.mbdyn.org/?Mailing_Lists) to the granted
access to the Gitlab repository while we work out a better solution

# MBDyn Developers Guidelines
MBDyn is developed primarily by _internal_ developers at Politecnico di 
Milano Department of Aerospace Science and Technology 
([DAER](http://www.aero.polimi.it/)).  
All contributions from _external_ developers however are welcome, of course.
As the project's maintainers, we only ask interested coders to follow the 
simple steps here presented.  
Since the main repository was moved to Git, we switched to the branching 
model described [here](https://nvie.com/posts/a-successful-git-branching-model/), 
so if you are a new developer, please read the page carefully prior to committing
and pushing changes.

## Guidelines for _external_ developers
As a general rule, before starting to write code, consider opening an 
[issue](https://gitlab.com/help/user/project/issues/index.md) to start a 
discussion with other developers. This, at the moment, **requires you to be
recognised as Project Member**, so please ask us to activate your account by 
writing to the [MBDyn users mailing list](https://www.mbdyn.org/?Mailing_Lists).

Developers that would like to contribute to MBDyn but are not in the project
regulars must:
 - fork the Gitlab repository hosted [here](https://gitlab.polimi.it/Pub/mbdyn.git)
 - checkout a fresh branch from the `develop` branch
 - commit and push to your branch in their forked repository
 - issue a merge request into the develop branch of MBDyn main repository

## Guidelines for _internal_ developers
Developers that would like to be included among the regulars in MBDyn must:

 - request a "developer" access permission to the MBDyn administrators, 
      @andomasarati in primis
 - clone the repository from [here](https://gitlab.polimi.it/Pub/mbdyn.git)
 - checkout a fresh branch from the `develop` branch
 - make sure that your contribution follows the code development model 
      indicated above before pushing
 - whenever they feel that their contribution would benefit from a discussion
      with other developers, issue a merge request instead of directly pushing
      to the main repository
    