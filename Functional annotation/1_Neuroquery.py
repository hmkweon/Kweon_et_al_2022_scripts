from neuroquery import fetch_neuroquery_model, NeuroQueryModel
from nilearn.plotting import view_img
import nibabel as nib
import numpy as np
from joblib import Parallel, delayed
import multiprocessing

LIST = np.loadtxt("./cognitive_atlas_concept.csv", dtype=str, delimiter=",") # csv file containing cognitive terms
encoder = NeuroQueryModel.from_data_dir(fetch_neuroquery_model())

def run_func(x):
  result = encoder(LIST[x])
  file = "./" + '_'.join(LIST[x].split()) + "_map"
  nib.save(result["brain_map"], file)      
    
Parallel(n_jobs=10)(delayed(run_func)(i) for i in range(0, LIST.size))

# some quries do not include every word in concpets. Need to drop them. 
# create individual word list
word_list = []
for i in range(0, LIST.size):
  word_list.extend(LIST[i].split())

word_list = np.unique(word_list)

for i in range(0, word_list.size):
  word_list.extend(LIST[i].split())


drop_list = []
for i in word_list:
  result = encoder(i)
  if result['raw_tfidf'].nnz == 0:
    drop_list.append(i)
    print(i)

# abductive
# chemonociception
# coreference
# dative
# diphthong
# down
# fixedness
# illocutionary
# indignation
# kinesthesia
# neologism
# numerosity
# of
# thermosensation
# to
# under

np.savetxt("./drop_list.csv", np.matrix(drop_list), fmt='%s')