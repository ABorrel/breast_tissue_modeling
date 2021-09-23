import pathFolder
import H295R_results

from os import path


### INIT FOLDERS
##################

PR_DATA = pathFolder.createFolder(path.abspath("../../") + "/DATA/")
PR_RESULTS = pathFolder.createFolder(path.abspath("../../") + "/RESULTS/")


c_H295 = H295R_results.H295R_results(PR_DATA, PR_RESULTS)
c_H295.main()
