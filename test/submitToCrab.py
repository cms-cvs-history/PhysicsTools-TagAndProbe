import os,sys
import string, re
from time import gmtime, localtime, strftime

## channels  = ["QCD_BCtoE_Pt20to30", "QCD_BCtoE_Pt30to80",
##              "QCD_BCtoE_Pt80to170", "QCD_EMenriched_Pt20to30",
##              "QCD_EMenriched_Pt30to80", "QCD_EMenriched_Pt80to170",
##              "TTbar", "Wenu", "Zee"]

## dataset  = ["/QCD_BCtoE_Pt20to30/Summer09-MC_31X_V3_SD_Ele15_QCD-v1/GEN-SIM-RECO",
##             "/QCD_BCtoE_Pt30to80/Summer09-MC_31X_V3_SD_Ele15-v1/GEN-SIM-RECO",
##             "/QCD_BCtoE_Pt80to170/Summer09-MC_31X_V3_SD_Ele15-v1/GEN-SIM-RECO",
##             "/QCD_EMEnriched_Pt20to30/Summer09-MC_31X_V3_SD_Ele15-v1/GEN-SIM-RECO",
##             "/QCD_EMEnriched_Pt30to80/Summer09-MC_31X_V3_SD_Ele15-v1/GEN-SIM-RECO",
##             "/QCD_EMEnriched_Pt80to170/Summer09-MC_31X_V3_SD_Ele15-v1/GEN-SIM-RECO",
##             "/TTbar/Summer09-MC_31X_V3_SD_Ele15-v1/GEN-SIM-RECO",
##             "/Wenu/Summer09-MC_31X_V3_SD_Ele15-v1/GEN-SIM-RECO",
##             "/Zee/Summer09-MC_31X_V3_SD_Ele15-v1/GEN-SIM-RECO"]


## dataset  = ["/QCD_BCtoE_Pt20to30/Summer09-MC_31X_V3_7TeV_SD_Ele15-v1/GEN-SIM-RECO",
##             "/QCD_BCtoE_Pt30to80/Summer09-MC_31X_V3_7TeV_SD_Ele15-v1/GEN-SIM-RECO",            
##             "/QCD_BCtoE_Pt80to170/Summer09-MC_31X_V3_7TeV_SD_Ele15-v1/GEN-SIM-RECO",            
##             "/QCD_EMEnriched_Pt20to30/Summer09-MC_31X_V3_7TeV_SD_Ele15-v1/GEN-SIM-RECO",            
##             "/QCD_EMEnriched_Pt30to80/Summer09-MC_31X_V3_7TeV_SD_Ele15-v1/GEN-SIM-RECO",
##             "/QCD_EMEnriched_Pt80to170/Summer09-MC_31X_V3_7TeV_SD_Ele15-v1/GEN-SIM-RECO",            
##             "/TTbar/Summer09-MC_31X_V3_7TeV_SD_Ele15-v1/GEN-SIM-RECO",            
##             "/Wenu/Summer09-MC_31X_V3_7TeV_SD_Ele15-v1/GEN-SIM-RECO",
##             "/Zee/Summer09-MC_31X_V3_7TeV_SD_Ele15-v1/GEN-SIM-RECO"]

## condor  = [1,0,0,1,1,0,0,1,1]
## condor  = [0,1,1,1,1,0,0,0,0]

## MyResilientArea = "/kalanand/Egamma_OctX"
## MyResilientArea = "/kalanand/Egamma_OctX_7TeV"

## MyConfigurationFile = "Electron_EDM_Ntuple_cfg.py"



## channels  = ["QCD_BCtoE_Pt20to30", "QCD_BCtoE_Pt30to80",
##              "QCD_BCtoE_Pt80to170", "QCD_EMenriched_Pt20to30",
##              "QCD_EMenriched_Pt30to80", "QCD_EMenriched_Pt80to170",
##              "Wenu", "Zee"]


## dataset  = ["/QCD_BCtoE_Pt20to30/e-gamma_ecal-GS_EGM1ELE15_v2-dc541c1fe0a4898b949d85e1df8ad863/USER",
##             "/QCD_BCtoE_Pt30to80/e-gamma_ecal-GS_EGM1ELE15_v2-dc541c1fe0a4898b949d85e1df8ad863/USER",            
##             "/QCD_BCtoE_Pt80to170/e-gamma_ecal-GS_EGM1ELE15_v2-dc541c1fe0a4898b949d85e1df8ad863/USER",            
##             "/QCD_EMEnriched_Pt20to30/e-gamma_ecal-GS_EGM1ELE15_v4-dc541c1fe0a4898b949d85e1df8ad863/USER",            
##             "/QCD_EMEnriched_Pt30to80/e-gamma_ecal-GS_EGM1ELE15_v4-dc541c1fe0a4898b949d85e1df8ad863/USER",
##             "/QCD_EMEnriched_Pt80to170/e-gamma_ecal-GS_EGM1ELE15_v3-dc541c1fe0a4898b949d85e1df8ad863/USER", 
##             "/Wenu/e-gamma_ecal-GS_EGM1ELE15_v3-dc541c1fe0a4898b949d85e1df8ad863/USER",
##             "/Zee/e-gamma_ecal-GS_EGM1ELE15_v2-dc541c1fe0a4898b949d85e1df8ad863/USER"]

## condor  = [0,0,0,0,0,0,0,0,0]

## MyResilientArea = "/kalanand/EwkWZ_OctX"
## MyConfigurationFile = "edm_ntpl_EwkWZ.py"
#prefix = "SD_Ele15_"


prefix = "Ntpl_"

channels  = ["Wenu"]
dataset  = ["/Wenu/Summer09-MC_31X_V3-v1/GEN-SIM-RECO"]

## channels  = ["Zee"]
## dataset  = ["/Zee/Summer09-MC_31X_V3-v1/GEN-SIM-RECO"]


condor  = [1]


            
MyResilientArea = "/kalanand/EwkWZ_OctX"
MyConfigurationFile = "ntpl_MC_test.py"



def changeMainConfigFile(outfile):
    fin  = open(MyConfigurationFile)
    pset_cfg      = "py_" + outfile + ".py"
    outfile_root  = prefix + outfile + ".root"
    fout = open(pset_cfg,"w")
    for line in fin.readlines():
        if  line.find("demo.root")!=-1:
            line=line.replace("demo.root",outfile_root)
        fout.write(line)
    print pset_cfg + " has been written.\n"


def changeCrabTemplateFile(outfile, index):
    fin  = open("crabTemplate.cfg")
    pset_cfg      = "py_" + outfile + ".py"
    pset_crab     = "crabjob_" + outfile + ".cfg"
    outfile_root  = prefix + outfile + ".root"
    fout = open(pset_crab,"w")
    for line in fin.readlines():
        if  line.find("mydataset")!=-1:
            line=line.replace("mydataset",dataset[index])
            #line += "\ndbs_url   = https://cmsdbsprod.cern.ch:8443/cms_dbs_ph_analysis_02_writer/servlet/DBSServlet\n"
        if line.find("myanalysis")!=-1:
            line=line.replace("myanalysis",pset_cfg)    
        if  line.find("myrootfile")!=-1:
            line=line.replace("myrootfile",outfile_root)
        if  line.find("myresilient")!=-1:
            line=line.replace("myresilient",MyResilientArea)    
        if line.find("glite")!=-1 and condor[index]!=0:
            line=line.replace("glite", "condor")        
        fout.write(line)        
    if condor[index]!=0:
        fout.write("ce_white_list = cmssrm.fnal.gov")

    print pset_crab + " has been written.\n"


    
###################
for i in range(len(channels)):
    changeMainConfigFile(channels[i])
    changeCrabTemplateFile(channels[i],i)

for i in range(len(channels)):
   submitcommand = "crab -create -cfg " + "crabjob_" + channels[i] + ".cfg"
   child   = os.system(submitcommand)
   child2   = os.system("crab -submit")
