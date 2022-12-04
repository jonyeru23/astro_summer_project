#The Second - Shussman

main paper : https://arxiv.org/abs/1610.05323

secondary paper: https://academic.oup.com/mnras/article/494/3/3927/5827327

Like the one before (i assume you first read SW16) the data is in apparent magnitude 
but the models gives us luminosity and temperature of the star at each time.

This time it is a bit more tricky, but you'll get it.

We don't assume a blackbody all the time, just for specific time (according to the secondary paper, after equation 5)

So we need to do two integrals, the first is to get the right value of L_nu in a specific filter
and the second is to correct it after Light Travel Time.

L_R_T.py concentrates most of the code for this specific model

get_model.py is the code for the model.

The parameters to fit are R500, M15, E51 and the offset as described in the papers. 

Because this is a bit more complex, there are extensive unittets for this model.

I encourage you to double test me and make sure all the graphs are by your standards.
