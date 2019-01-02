# 02450_DeepLearningProject
Deep Learning Project for the class 02450 DeepLearning taught at Denmarks Technical University (DTU) Fall 2018

For repoucing the results the following has to be done:

1.) Prepare the databases for the Training input by running the script <b>merge_db.py</b>. The <b>mixdim.db</b> containing all the structure files can be downloaded from the CMR repository at https://cmr.fysik.dtu.dk/mixdim/mixdim.html#mixdim. The python script will pick all Crystalography Open Database entries (http://www.crystallography.net/cod/) and assign the correct class. The definition of the class is explained in the underlying paper: Haastrup, Sten, et al. “Definition of a scoring parameter to identify low-dimensional materials components” arXiv preprint arXiv:1806.03173 (2018).

2.) After the databases are prepared, the fingerprints need to be created by running the <b>prepare_fingerprints.py</b> script.

3.) The main results can then be reproduced using the ClassifyUseCNN notebook. It also includes a more detailed explanation on the classification scheme.
