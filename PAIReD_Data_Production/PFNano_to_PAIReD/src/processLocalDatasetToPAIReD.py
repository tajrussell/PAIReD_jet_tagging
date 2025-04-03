
"""make N-tuples of PAIReD jets from MC data set

This script processes an entire data set from DAS to the PAIReD
data format.
"""
from dbs.apis.dbsClient import DbsApi
from processFileToPAIReD import makeNtuplesPAIReDjointMC
import argparse
import multiprocessing as mp
import functools
import os.path
import os


def process(inputpath, outputpath, batchsize, N_update,
            N_events_to_process, physics_process, PAIReD_geometry, inputfile=None):
    """
    Processes the file with index i.

    Args:
        i (int): index of the file to be processed
    """
    # path to PFNanoAOD file
    PFNanopath = inputpath + inputfile
    # get index of the PFNanoAOD file and define the output file path
    i = inputfile.rpartition("_")[-1][:-5]
    PAIReDpath = outputpath + "/PAIReD_%s.root" % i
    # check if file even exists
    try:
        makeNtuplesPAIReDjointMC(PFNanopath, PAIReDpath, batchsize = batchsize,
                                N_update = N_update, N_events_to_process = N_events_to_process, 
                                physics_process=physics_process, PAIReD_geometry=PAIReD_geometry)
    except FileNotFoundError:
        print("\n** file", i, "does not exist. Therefore, skip it!\n")
        return 0


def processDataset(dataset,
                   outputpath,
                   inputpath,
                   batchsize,
                   first_file,
                   N_files,
                   N_update,
                   N_events_to_process,
                   N_threads,
                   physics_process,
                   PAIReD_geometry):

    print("\n#################################################################")
    print("  Start processDataset()")
    print("#################################################################\n")

    print("** Process the following data set:")
    print("   ", dataset)
    print("** input path:", inputpath)

    filelist = [dataset+f for f in os.listdir(inputpath+dataset) if os.path.isfile(os.path.join(inputpath+dataset, f))]

    print("** Number of files in the data set:", len(filelist))
    if N_files==-1:
        N_files = len(filelist) - first_file
    last_file = min(first_file + N_files, len(filelist))
    print("** Number of files to be processed:", last_file-first_file)
    print("** Start at file:", first_file)

    print("\nSettings for processFileToPAIReD:")
    print("** Batch size:", batchsize)
    print("** N_update:", N_update)
    print("** Number of events to process:", N_events_to_process)
    print("** PAIReD geometry:", PAIReD_geometry)
    print("** Create output directory if it does not exist yet:", outputpath)
    os.makedirs(outputpath, exist_ok=True)

    # if not run in parallel
    if N_threads == 1:
        print("** Use only one thread")
        print("\n** Start processing the data set\n\n")
        for inputfile in filelist[first_file:last_file]:
            process(inputpath, outputpath, batchsize, N_update,
                        N_events_to_process, physics_process, PAIReD_geometry, 
                        inputfile=inputfile)
    # if run in parallel
    else:
        process_for_pool = functools.partial(process, inputpath,
                                            outputpath, batchsize, N_update,
                                            N_events_to_process, physics_process, PAIReD_geometry)
        print("** Use %i threads in parallel" % N_threads)
        print("\n** Start processing the data set\n\n")
        with mp.Pool(processes=N_threads) as pool:
            pool.map(process_for_pool, filelist[first_file:last_file])

    print("\n\n** Done processing\n")

    print("\n#################################################################")
    print("  End processDataset()")
    print("#################################################################\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process MC PFNano data set to the PAIReD data format for training.")
    parser.add_argument("dataset", help="Name of the data set (NanoAOD processed by PFNano, has to contain particle flow candidates (PFcands))")
    parser.add_argument("outputpath",  help="Path for the ROOT output files")
    parser.add_argument("-i", "--inputpath", default="/home/trussel1/bruxhcc/", help="path to where the dataset is held locally")
    parser.add_argument("-b", "--batchsize", default="100 MB", help="Size of the batches of events that are handled simultaneously (default is '100 MB'). Can be an integer (=number of events) or a string (specifying the data size, e.g. '100 MB')")
    parser.add_argument("-u", "--N-update", type=int, default=10, help="Number of batches processed before printing an update in the terminal: e.g., 'process event batch 10 of 15' (default is 10)")
    parser.add_argument("-n", "--N-files", type=int, default=-1,  help="Number of files to be processed. -1 means all (default is -1)")
    parser.add_argument("--start-file", type=int, default=1,  help="Index of first file to be processed (default is 0)")
    parser.add_argument("-e", "--N-events", type=int, default=-1,  help="Number of events to be processed per file. -1 means all (default is -1)")
    parser.add_argument("-t", "--N-threads", type=int, default=1,  help="Number of threads (default is 1)")
    parser.add_argument("--physicsprocess", type=int, default=0,  help="Integer indicating the physics process. Default is 0.")
    parser.add_argument("-g", "--geometry", type=str, default="Ellipse",  help="String indicating the geometry of the PAIReD jet. Default is Ellipse.")
    args = parser.parse_args()

    processDataset(args.dataset, args.outputpath,
                   inputpath=args.inputpath,
                   batchsize=args.batchsize,
                   first_file=args.start_file,
                   N_files=args.N_files,
                   N_update=args.N_update,
                   N_events_to_process=args.N_events,
                   N_threads=args.N_threads,
                   physics_process=args.physicsprocess,
                   PAIReD_geometry=args.geometry)
