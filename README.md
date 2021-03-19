# Assignment3
This is the repo for assignment 3. It is a two part assignment and is due at 1pm on April 5th.

## Part 1 - Sampling Gamma random variable

Use the IndeptGamma.Rmd script to use MCMC to produce samples from a Gamma(2.3,2.7) random variable.

Show how the performance deteriorates when you run the algorithm to produce samples from a Gamma(0.1,0.01).

Use the Gelman plot to illustrate the deterioration in performance.

Find a proposal kernel (i.e. q(x->x’), the way of producing new candidate values) that performs more efficiently.
Again use the Gelman plots to show how performance has changed. Discuss how the performance has improved



## Part 2 - Code Breaking

Use one of the methods you have met in the course to try to break the coded message one of the following files:

* CodedMessage_Short.txt
* CodedMessage_Med.txt

They are two extracts from a famous novel, the Med one is longer than the Short one. The longer message will give you a beter chance of breaking the code, but will run more slowly. (I suggest starting with the short one, and once you are happy with your code trying it on the med one.)

It may not be easy to break the code completely without resorting to manual ‘tweaks’ at the end, but you should be able to get to the point at which you can work out what the text is saying.

### Clues

You are given the following clues:

* It is a simple substitution code
e.g. “a/A” might stand for “m/M”;  “b/B” might stand for “f/F”
* The substitution is the same for the entire message.
* You are given a reference table of pairwise letter frequencie extracted from 7 novels (see the file LetterPairFreqFrom7Novels.txt)
* We are doing this just after we learned about optimization and MCMC.


Write both parts up as an Rmarkdown file and include a description of the methods you are using (for both parts of the assigment) and a discussion of your results.
Due April 5th at 1pm. Please commit them to your own versions of this repository (i.e., the ones for which you received an invite to create from Github classroom.)
