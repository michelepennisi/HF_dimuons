Macro to lauch powheg simulations using pythia6 

The aim of these macros is to perform full ALICE MC simulation either on the local machine or in the grid. Note that for both modes, an alien connection (with alien-token-init) must be established.

The main files are:
-simrun.sh: This is the main script to be run in the case of local testing. It then runs respectively the chain of other macros needed for the simulation (generation + reconstruction + filtering (ESD to AOD).

***************************************************************************
To run simrun.sh locally:
./simrun.sh --run <runNumber> --eventsPerJob <number of event to be generated> --energy <collision energy> --transport <geant3 or geant4> --generator <pythia or powheg> --collisionSystem <pp or nn or pn or np> --beamConfig <pPb or Pbp>

Simrun.sh will call the following macros:
-Config.C: this macro will set the configuration needed for the simulation. It calls all the aliroot library needed. It will set the energy collision and start the primary particle generator (Pythia or Powheg). It also calls the transport code (geant3 or geant4).

-sim.C: this macro will start the particle generation and it is also responsible for setting the right paths for the OCDB and the detector needed files (like the alignment). After the execution of this macro, the kinematics of the MC particles are stored in a root file called (Kinematics.root) and the detector hits are stored in ?hits.root. Those files will be used by the reconstruction macro (rec.C)

-rec.C: is responsible for the reconstruction. The good alignment path for the studied period must be given. The result of this macro is an ESD file "AliESDs.root".

-AODtrainsim.C: This macro will finally filter the ESDs into AODs. The result of this macro is "AliAOD.Muons.root

***************************************************************************

To run the simulation on the grid:
Three other files must be used here.
-Submit.C: this will ensure the connection to alien and submit the jobs to the grid.
-run.jdl: This file contains the settings for the grid analysis (aliroot version, time price of the job?.) and the alien path to the files needed to run the simulation (all the files used for the local simulation).
-runList.txt: Submit.C will submit the jobs for the run that are included in runList.txt. The structure of this file must not be changed. Each line contains a run number, the number of jobs per run, the number of events per job, the transport software (geant3 or geant4), the particle generator (pythia or powheg),  the collisions system (pp or pn or np or nn) and the beam configuration (pPb or Pbp).
For example the line "266304 2 10 geant3 powheg pp pPb" will submit 2*10 events for the run 266304, simulating a p-Pb collision for p-p binary collisions, simulated with powheg and transported with geant3.

This script is not fully automatic so you need to make sure that you copied the files mentioned under "InputFile" in the run.jdl to your alien working directory. Also make that all the paths mentioned in run.jdl must corresponds to your own working directory and containing your cern username
***************************************************************************

*The files: g4libs.C and basiclibs.c are needed for geant4 and you do not need to modify them.
*the file base_powheg.input contains the settings for the powheg generator.
