import os

def main():
    os.system("cp ../AliAnalysisTaskDimuon_HighMass.cxx .")
    os.system("cp ../AliAnalysisTaskDimuon_HighMass.h .")
    os.system("cp ../AddTaskDimuon_HighMass_Local.C .")
    
    os.system("aliroot -q ../ReadMCDimuon_HighMass_Local.C+")

main()