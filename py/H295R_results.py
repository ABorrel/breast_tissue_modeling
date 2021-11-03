from os import path
import toolbox
import pathFolder
import runExternal

class H295R_results:
    """
    Class use to load and parse table from Kaumaus et al 2006, study develop from H295R
    """

    def __init__(self, pr_data, pr_results):

        pr_Karmaus2016 = pr_data + "karmaus_2016_SI/kfw002_Supplementary_Data/"
        pr_Haggard2018 = pr_data + "haggard_2018_SI/Haggard_et_al_Supplemental_Files_Revision1_30Oct2017/"

        # extract data from cardona 2021
        p_process_data = pr_data + "e2up_p4up_Cardona_2021_SI.xlsx"
        self.p_process_data = p_process_data
        
        # 1- load results
        self.d_E2up = toolbox.loadExcelSheet(self.p_process_data, name_sheet='H295R_E2up', k_head = "CASN_Protect")
        self.d_P4up = toolbox.loadExcelSheet(self.p_process_data, name_sheet='H295R_P4up', k_head = "CASN_Protect")


        pr_results = pathFolder.createFolder(pr_results + "H295R/")
        self.pr_results = pr_results

        # load from kaumas
        self.pf_chem_MTT = pr_Karmaus2016 + "toxsci-15-0570-File006.csv"
        self.pf_concentration_response = pr_Karmaus2016 + "toxsci-15-0570-File006.csv"
        self.pf_hormone_resp = pr_Karmaus2016 + "toxsci-15-0570-File008.csv"
        self.pf_hormone_hitc = pr_Karmaus2016 + "toxsci-15-0570-File009.csv"
        self.pf_hormone_CR = pr_Karmaus2016 + "toxsci-15-0570-File010.csv"
        self.pf_raw_hormone_data = pr_Karmaus2016 + "toxsci-15-0570-File011.csv"
        self.pf_hormone_change = pr_Karmaus2016 + "toxsci-15-0570-File014.csv"
        
        # load raw data from haggard
        self.pf_raw_data = pr_Haggard2018 + "Supp3_H295R_master_table_2017-08-08.csv"
        self.pf_raw_anova = pr_Haggard2018 + "Supp4_OECD_GLOBAL_ANOVA_output_pValues2017-08-09.txt"

        # define hormones for the analysis
        self.l_hormones = ["OHPREG", "PROG", "OHPROG", "DOC", "CORTISOL", "X11DCORT", "ANDR", "TESTO", "ESTRONE", "ESTRADIOL"]

    def loadsinglehitc(self):
        """With Kaumas data"""

        l_lines_hitc = toolbox.loadMatrixToList(self.pf_hormone_hitc, sep = ",")
        d_hitc = {}
        if not "d_chem" in self.__dict__:
            self.loadChemicalsMapping()
        
        for line_hitc in l_lines_hitc:
            chid = line_hitc["chid"]
            if chid == "NA":
                continue
            aenm = line_hitc["aenm"]
            hormone = aenm.split("_")[2]
            up_down = aenm.split("_")[-1]
            if not chid in list(d_hitc.keys()):
                d_hitc[chid] = {}
            if line_hitc["hitc"] != "0":
                if up_down == "up":
                    d_hitc[chid][aenm] = float(line_hitc["max_med"])
                else:
                    d_hitc[chid][aenm] = -float(line_hitc["max_med"])
        
        self.d_hitc = d_hitc

    def summaryHitc(self):
        """With Kaumas data"""

        p_filout = self.pr_results + "single_hitc_matrix.csv"
        if path.exists(p_filout):
            d_out = toolbox.loadMatrix(p_filout)
            self.d_single_hit = d_out
            return 
        
        if not "d_hitc" in self.__dict__:
            self.loadsinglehitc()
            
        d_out = {}
        for chem in self.d_hitc.keys():
            if not chem in list(d_out.keys()):
                d_out[chem] = {}
                for h in self.l_hormones:
                    d_out[chem][h] = 0
            
            for aenm in self.d_hitc[chem].keys():
                hormone = aenm.split("_")[2]
                hitc = self.d_hitc[chem][aenm]
                if hitc != 0:
                    d_out[chem][hormone] = hitc
        
        filout = open(p_filout, "w")
        filout.write("CASRN\t%s\tSUM_HORMONE_CHANGED\n"%("\t".join(self.l_hormones)))
        for chem in d_out.keys():
            casrn = self.d_chem_mapping[chem]["casn"]
            nb_h_changed = 0
            for h in self.l_hormones: 
                if d_out[chem][h] != 0: 
                    nb_h_changed = nb_h_changed +1 
            filout.write("%s\t%s\t%s\n"%(casrn, "\t".join([str(d_out[chem][h]) for h in self.l_hormones]), nb_h_changed))
        filout.close()
    
        # write a summary
        runExternal.barplotHormones(p_filout)

        d_out = toolbox.loadMatrix(p_filout)
        self.d_single_hit = d_out

    def loadChemicalsMapping(self):
        """With Kaumas data"""

        d_chem = {}
        l_lines = toolbox.loadMatrixToList(self.pf_hormone_resp, sep = ",")

        for line_file in l_lines:
            chid = line_file["chid"]
            spid = line_file["spid"]
            casn = line_file["casn"]
            chnm = line_file["chnm"]

            if not chid in list(d_chem.keys()) and chid != "NA":
                d_chem[chid] = {}
                d_chem[chid]["chid"] = chid
                d_chem[chid]["spid"] = spid
                d_chem[chid]["casn"] = casn
                d_chem[chid]["chnm"] = chnm
        
        self.d_chem_mapping = d_chem

    def loadChemicalCR(self):
        """With Kaumas data"""


        p_filout = self.pr_results + "CR_hitc_matrix.csv"
        if path.exists(p_filout):
            self.d_CR_hit = toolbox.loadMatrix(p_filout)
            return
        
        if not "d_chem_mapping" in self.__dict__:
            self.loadChemicalsMapping

        d_out = {}
        
        l_lines_CR = toolbox.loadMatrixToList(self.pf_hormone_change, sep = ',')
        for line_CR in l_lines_CR:
            chid = line_CR["chid"]
            casn = self.d_chem_mapping[chid]["casn"]
            d_out[casn] = {}

            for h in self.l_hormones:
                d_out[casn][h] = line_CR[h]
        
        filout = open(p_filout, "w")
        filout.write("CASRN\t%s\tSUM_HORMONE_CHANGED\n"%("\t".join(self.l_hormones)))
        for casrn in d_out.keys():
            nb_h_changed = 0
            for h in self.l_hormones: 
                if d_out[casrn][h] != "0": 
                    nb_h_changed = nb_h_changed +1 
            filout.write("%s\t%s\t%s\n"%(casrn, "\t".join([str(d_out[casrn][h]) for h in self.l_hormones]), nb_h_changed))
        filout.close()
        runExternal.barplotHormones(p_filout)
        self.d_CR_hit = toolbox.loadMatrix(p_filout)

    def corHormoneEndpoint(self):
        """With Kaumas data"""


        pr_out = pathFolder.createFolder(self.pr_results + "CorHorm/")
        p_filout = pr_out + "horm_resp.csv"
        filout = open(p_filout, "w")
        filout.write("CASRN\t%s\n"%("\t".join(self.l_hormones)))


        # plot correlation for 2 hormones
        for chem in self.d_CR_hit.keys():
            filout.write("%s\t%s\n"%(chem, "\t".join([self.d_CR_hit[chem][h] for h in self.l_hormones])))
        filout.close()     

        # plot correlation
        runExternal.corHormResponse(p_filout, pr_out)

    def drawResponseCurve(self, l_hormones):

        pr_out = pathFolder.createFolder(self.pr_results + "Responce_curve/")
        l_d_raw = toolbox.loadMatrixToList(self.pf_raw_data, sep = "\t")

        l_d_anova = toolbox.loadMatrixToList(self.pf_raw_anova, sep = "\t")

        d_out = {}
        for hormone in l_hormones:
            d_out[hormone] = {}
        for d_raw in l_d_raw:
            casrn  = d_raw["casn"]
            chnm = d_raw["chnm"]
            hormone = d_raw["steroid"]
            
            if hormone in l_hormones:
                if not casrn in list(d_out[hormone].keys()):
                    d_out[hormone][casrn] = {}
                    d_out[hormone][casrn]["name"] = chnm
                    d_out[hormone][casrn]["concentration"] = []
                    d_out[hormone][casrn]["response"] = []
                    d_out[hormone][casrn]["anova"] = []
                
                d_out[hormone][casrn]["concentration"].append(d_raw["conc"])
                d_out[hormone][casrn]["response"].append(d_raw["uM"])

                for d_anova in l_d_anova:
                    if d_anova["casn"] == casrn:
                        #print("++++++++++++++++")
                        #print(d_raw["conc"])
                        #print(d_anova["conc"])
                        if round(float(d_anova["conc"]), 0) == round(float(d_raw["conc"]), 0):
                             d_out[hormone][casrn]["anova"].append(d_anova[hormone])
        
        for hormone in d_out.keys():
            pr_byhormone = pathFolder.createFolder(pr_out + hormone + "/")
            for chem in d_out[hormone].keys():
                if chem == "":
                    continue
                if hormone == "Estradiol" and chem in list(self.d_E2up.keys()):
                    eff_pot = self.d_E2up[chem]["Efficacy/potency"]
                elif hormone == "Progesterone" and chem in list(self.d_P4up.keys()):
                    eff_pot = self.d_P4up[chem]["Efficacy/potency"]
                else:
                    eff_pot = "NA"
                
                try:fold_change = self.d_CR_hit[chem][hormone.upper()]
                except: fold_change = 0

                pr_eff_pot = pathFolder.createFolder(pr_byhormone + eff_pot.replace(" ", "_") + "/")
                p_filout = pr_eff_pot + chem
                filout = open(p_filout, "w")
                filout.write("Concentration\tResponse\tchnm\tCR results\tEfficacy/potency\tAnova\n")

                i = 0
                imax = len(d_out[hormone][chem]["concentration"])
                while i < imax:
                    filout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(d_out[hormone][chem]["concentration"][i], d_out[hormone][chem]["response"][i], d_out[hormone][chem]["name"], fold_change, eff_pot, d_out[hormone][chem]["anova"][i]))
                    i = i + 1
                filout.close()

                runExternal.drawResponseCurve(p_filout, hormone, 0.05)





    def main(self):
        """
        Load data from stereogenesis 
        """

        self.loadChemicalsMapping()
        self.loadChemicalCR()
        self.drawResponseCurve(["Estradiol"])
        #self.loadsinglehitc()
        #self.summaryHitc()
        #self.loadChemicalCR()

        #self.corHormoneEndpoint()

