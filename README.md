# astro_summer_project
In this repository, I am analyzing the data from supernova SN2016gkg.

I am combining data from two groups of researches, and it is centrelized in the excel files in this repository.

The main objective is to identify at what stage of the supernova's explosion we started recieving data and disproving (maybe) the following paper, 
claming the data is from Shock breakout stage: https://www.nature.com/articles/nature25151

I am using two main papers to do so.

The first is based on Iair Arcavi's paper : https://ui.adsabs.harvard.edu/abs/2017ApJ...837L...2A/abstract
Which analizes the same supernova whithout the data from the other group.

Arcavi used the model described in this paper : https://arxiv.org/abs/1607.03700 

For the second model I am using is these two papers: 

https://arxiv.org/abs/1610.05323
which gives the main aspects of the model

and this one, which corrects the first one regarding light travel time: https://academic.oup.com/mnras/article/494/3/3927/5827327

The main libraries i am using and a worth a read: Astropy, Pysynphot, emcee, nnumpy and pandas.

further explanations are in the Readme.md files in the folders and the files themselfs.

I wish you luck! and hope you get to the true answer regarding this supernova!

Important - using pysynphot is a real pain in the ass. 

you need to download a lot of files just like they say and set an ENVIRONMENT variable
accordingly. When you run it on the nove, you need to set the environment variable every time you run it.
My recomendation is to download all the files in a place you can find them easily and put the command to set the environment
variable there as well. Link to installing pysynphot: https://pysynphot.readthedocs.io/en/latest/