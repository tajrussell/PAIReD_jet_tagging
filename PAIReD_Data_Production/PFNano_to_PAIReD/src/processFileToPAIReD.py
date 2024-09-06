"""make N-tuples of PAIReD jets from MC data file

This script allows the user to convert MC simulations in the format of 
PFNano files to ROOT files that can be used in the PAIReD tagger network 
for training.

This tool accepts ROOT files (.root) processed by PFNano.

This script requires the following modules to be installed within the
Python environment you are running this script in:

    * uproot
    * awkward
    * numpy
    * tools.helpers

This file can also be imported as a module and contains the following
functions:

    * makeNtuplesPAIReDjointMC [main] - converts the data input file
"""


# import all modules needed
import numpy as np
import awkward as ak
import uproot
import vector
vector.register_awkward()
from tools.processEvents import processEvents
from tools.branchnames import outputTreeType
import json
import time
import argparse
import os

# ******************************************************************************


"""
 * Function: makeNtuplesPAIReDjointMC  [main function]
 * -----------------------------------------------------
Performs the data translation from PFNano output files (extended NanoAODs)
to a data format that the PAIReD tagger network can read in and train with.
Remark: this function is designed to deal with CMS Monte Carlo simulations and
NOT real CMS data.

Parameters
----------
inputFilePath : str
    Path to the input file (NanoAOD processed by PFNano, has to contain particle
    flow candidates (PFcands))
outputFilePath : str
    Path and name of the output file
batchsize : int, optional (default is "100 MB")
    Size of the batches of events that are handled simultaneously. 
    If int = number of events,
    if str = data size, e.g., "100 MB"
N_update : int, optional (default is 10)
    Number of batches to process before printing an update of the form:
        process event batch 10
N_events_to_process : int, optional (default is -1)
    Number of events to process. =-1 means all events are processed.
physics_process : int, optional (default is 0)
    Integer indicating the physics process in the file. =0 means unknown/not specified.
PAIReD_geometry : str, optional (default is "Ellipse")
    String indicating the geometry of the PAIReD jet.


Returns
-------

    nothing, but creates an output file with the transformed data
"""

def makeNtuplesPAIReDjointMC(inputFilePath, outputFilePath, batchsize = "100 MB",
    N_update = 10, N_events_to_process = -1, physics_process = 0, PAIReD_geometry = "Ellipse"):

    print("\n*****************************************************************")
    print("  Start makeNtuplesPAIReDjointMC()")
    print("*****************************************************************\n")

    print("** Open input file:", inputFilePath)
    # Access the "Events" tree of the input file
    with uproot.open(inputFilePath) as inputFile:
        inputTree = inputFile['Events']

        print("** Create output file:", outputFilePath)
        with uproot.recreate(outputFilePath,
                             compression=uproot.LZ4(4)) as outputFile:


            # create ROOT Tree in output file
            outputTree = outputFile.mktree("tree", outputTreeType, 
                                           title="PAIReD jets for training")

            print('** Create tree "tree" entitled "PAIReD jets for training"',
                    'in output file')

            # number of events in input file
            N_events = inputTree.num_entries
            print("** Total number of events:", N_events)
            
            if N_events_to_process > N_events or N_events_to_process==-1:
                N_events_to_process = N_events
                print("** Process all events")
            else: 
                print("** Process first %i events only" % N_events_to_process)
            
            print("** Batch size:", batchsize)
            print("** Physics process:", physics_process)
            print("** PAIReD geometry:", PAIReD_geometry)

            # names of the branches that will be used from the input tree
            script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
            rel_path = 'tools/input_branches.json'
            json_file_path = os.path.join(script_dir, rel_path)
            with open(json_file_path) as json_file:
                inputBranchNames = json.load(json_file)
            
            
            t_start = time.time()
            print("\n** Start processing data")

            # number of processed batches
            N_processed = 0

            # iterate over all batches
            for Events in inputTree.iterate(expressions=inputBranchNames,
                    cut="(nJet>=2) & (nGenJet>=1)", step_size=batchsize, entry_stop=N_events_to_process):
                
                # output update every N_update batches
                if N_processed % N_update == 0:
                    print("process event batch", N_processed)
                
                # veto if there are no events with enough jets to create a
                # jet pair in the batch
                if len(Events) == 0:
                    print("   !!! No events with more than one jet in batch !!!")
                    print("   * Batch %i was not processed" % N_processed)
                    N_processed += 1
                    continue

                # process the events
                DataPAIReD = processEvents(Events, physics_process=physics_process, PAIReD_geometry=PAIReD_geometry)
                # check if processing worked
                if DataPAIReD == False:
                    print("   * Batch %i was not processed" % N_processed)
                    N_processed += 1
                    continue
                
                # fill in the tree
                outputTree.extend(DataPAIReD)

                N_processed += 1

            t_end = time.time()
            dt = t_end - t_start

            print("\n** Written", outputTree.num_entries,
                    "PAIReD jets to the output file in %.2f seconds" % dt)
            print("** Computing time per PAIReD jet:",
                    "%.3f milliseconds" % (dt / outputTree.num_entries *1e3))
            print("** Computing time per event:",
                    "%.3f milliseconds" % (dt / N_events_to_process *1e3))

    print("\n*****************************************************************")
    print("  End makeNtuplesPAIReDjointMC()")
    print("*****************************************************************\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process MC PFNano files to the PAIReD data format for training.")
    parser.add_argument("inputFilePath", help="Path to the ROOT input file (NanoAOD processed by PFNano, has to contain particle flow candidates (PFcands))")
    parser.add_argument("outputFilePath",  help="Path and name of the ROOT output file")
    parser.add_argument("--batchsize", default="100 MB", help="Size of the batches of events that are handled simultaneously (default is '100 MB'). Can be an integer (=number of events) or a string (specifying the data size, e.g. '100 MB')")
    parser.add_argument("--N-update", type=int, default=10, help="Number of batches processed before printing an update in the terminal: e.g., 'process event batch 10 of 15' (default is 10)")
    parser.add_argument("-n", "--nevents", type=int, default=-1,  help="Number of events to be processed")
    parser.add_argument("-p", "--physicsprocess", type=int, default=0,  help="Integer indicating the physics process. Default is 0.")
    parser.add_argument("-g", "--geometry", type=str, default="Ellipse",  help="String indicating the geometry of the PAIReD jet. Default is Ellipse.")
    args = parser.parse_args()

    makeNtuplesPAIReDjointMC(args.inputFilePath, args.outputFilePath,
        batchsize=args.batchsize, N_update=args.N_update, 
        N_events_to_process=args.nevents, physics_process=args.physicsprocess, PAIReD_geometry=args.geometry)
