from os import system, path, remove, chdir, getcwd, listdir, name
import subprocess 


P_RSCRIPTS = "../R/"
PR_BIOTRANSFORMER = "C:/Users/Aborrel/research/Silent_Spring/PFAS/BioTransformerJar/biotransformerjar"

R_BIN = "C:\\Program Files\\R\\R-4.0.2\\bin\\Rscript.exe"
P_RQSAR_linux = "/mnt/c/Users/AlexandreBorrel/research/development/QSAR-QSPR"
P_RQSAR_window = "c:/Users/aborr/research/development/QSAR-QSPR/"
######
# Main functions

def runRCMD(cmd, out = 0):

    workdir = getcwd()
    chdir(P_RSCRIPTS)
    if name == "nt":
        l_elem = cmd.split(" ")
        cmd_line = [R_BIN] + l_elem
        print(cmd_line)
        p = subprocess.Popen(cmd_line)
        (output, err) = p.communicate() 
        p.wait()
        print(err)
    else:
        print(cmd)
        system(cmd)
    chdir(workdir)



# for the Kaumas dataset

def corHormResponse(p_filin):
    cmd = "./corHormResponse.R %s"%(p_filin)
    runRCMD(cmd)

def barplotHormones(p_filin):
    cmd = "./barplotHormones.R %s"%(p_filin)
    runRCMD(cmd)   

def drawResponseCurve(p_filin, hormone):
    cmd = "./draw_RC.R %s %s"%(p_filin, hormone)
    runRCMD(cmd)  