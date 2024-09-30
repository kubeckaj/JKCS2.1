#TODO:

# We want to start using https://xacs.xmu.edu.cn/docs/mlatom/index.html

# those are two functions for training and evaluating that I at least want to have here
#   def training_nn
#   def evaluating_nn 
# find inspiration in SchNetPack.py

#there must be many modifications to arguments.py and also a bit to ../JKML.py

#I will forward you structures from JKML.py and they will be in ASE format so you will need to convert them first
#I would use e.g. classmethodfrom_numpy(coordinates: ndarray, species: ndarray)→ molecule

# I will be at the same time working on QML.py to reduce the weight of ./JKML.py

#there will likely be necessary many Git pulls/pushes
#I do suggest that we learn a bit more about GitHub, create another branch where you will have full access and we will do the modifications
#when all done, I merge the branches and push to the main 

#there might oc some obstacle with installation too, I can help with that

#I want to add nearly all models. I do not care about global descriptiors for now. PhysNet and ANI should be there: https://xacs.xmu.edu.cn/docs/mlatom/tutorial_mlp.html#tutorial-mlp-supported

#later, I will want to have access to all pretrained models: https://xacs.xmu.edu.cn/docs/mlatom/tutorial_universalml.html
