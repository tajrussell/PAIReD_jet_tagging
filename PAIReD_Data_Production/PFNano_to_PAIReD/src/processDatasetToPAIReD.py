"""make N-tuples of PAIReD jets from MC data set

This script processes an entire data set from DAS to the PAIReD
data format.
"""
from processFileToPAIReD import makeNtuplesPAIReDjointMC
import argparse
import multiprocessing as mp
import functools
import os.path
import os

def process(storageprefix, datasetpath, filename, outputpath, batchsize, N_update,
            N_events_to_process, physics_process, PAIReD_geometry, i=0):
    """
    Processes the file with index i.

    Args:
        i (int): index of the file to be processed
    """
    # catch files named PFNANOAOD*_1-j.root by using indices i=1000+j instead
    if i<1000:
        PFNanopath = storageprefix + datasetpath + "/" + filename + "%i.root" % i
    else:
        PFNanopath = storageprefix + datasetpath + "/" + filename + "1-%i.root" % (i-1000)
    PAIReDpath = outputpath + "/PAIReD_%i.root" % i

    # check if file even exists
    try:
        makeNtuplesPAIReDjointMC(PFNanopath, PAIReDpath, batchsize = batchsize,
                                N_update = N_update, N_events_to_process = N_events_to_process, 
                                physics_process=physics_process, PAIReD_geometry=PAIReD_geometry)
    except FileNotFoundError:
        print("\n** file", i, "does not exist. Therefore, skip it!\n")
        return 0


def processDataset(datasetpath,
                   outputpath,
                   storageprefix = "root://grid-cms-xrootd.physik.rwth-aachen.de/",
                   filename = "PFNANOAODSIM_allPF_allBTV_",
                   batchsize = "100 MB",
                   first_file = 1,
                   N_files = 1,
                   N_update = 10,
                   N_events_to_process = -1,
                   N_threads = 1,
                   physics_process = 0,
                   PAIReD_geometry = "Ellipse"):

    print("\n#################################################################")
    print("  Start processDataset()")
    print("#################################################################\n")

    print("** Process data set stored in the following directory:")
    print("   ", storageprefix + datasetpath)
    print("** Number of files to be processed:", N_files)

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
        for i in range(first_file, first_file + N_files):
            process(storageprefix, datasetpath, filename, outputpath, batchsize, N_update,
                        N_events_to_process, physics_process, PAIReD_geometry, i=i)
    # if run in parallel
    else:
        process_for_pool = functools.partial(process, storageprefix, datasetpath,
                                            filename, outputpath, batchsize, N_update,
                                            N_events_to_process, physics_process, PAIReD_geometry)
        print("** Use %i threads in parallel" % N_threads)
        print("\n** Start processing the data set\n\n")
        with mp.Pool(processes=N_threads) as pool:
            pool.map(process_for_pool, range(first_file, first_file + N_files))

    print("\n\n** Done processing\n")

    print("\n#################################################################")
    print("  End processDataset()")
    print("#################################################################\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process MC PFNano data set to the PAIReD data format for training.")
    parser.add_argument("dataset", help="Path of the files of the data set to process (NanoAOD processed by PFNano, has to contain particle flow candidates (PFcands))")
    parser.add_argument("outputpath",  help="Path for the ROOT output files")
    parser.add_argument("-p", "--storageprefix", default="root://grid-cms-xrootd.physik.rwth-aachen.de/", help="Prefix of the storgae site the data set is stored on (default for T2_DE_RWTH site)")
    parser.add_argument("-f", "--filename", default="PFNANOAODSIM_allPF_allBTV_", help="Default filename for files in data set without index number and root extension at the end (default is 'PFNANOAODSIM_allPF_allBTV_')")
    parser.add_argument("-b", "--batchsize", default="100 MB", help="Size of the batches of events that are handled simultaneously (default is '100 MB'). Can be an integer (=number of events) or a string (specifying the data size, e.g. '100 MB')")
    parser.add_argument("-u", "--N-update", type=int, default=10, help="Number of batches processed before printing an update in the terminal: e.g., 'process event batch 10 of 15' (default is 10)")
    parser.add_argument("-n", "--N-files", type=int, default=1,  help="Number of files to be processed")
    parser.add_argument("--start-file", type=int, default=1,  help="Index of first file to be processed (default is 1)")
    parser.add_argument("-e", "--N-events", type=int, default=-1,  help="Number of events to be processed per file")
    parser.add_argument("-t", "--N-threads", type=int, default=1,  help="Number of threads (default is 1)")
    parser.add_argument("--physicsprocess", type=int, default=0,  help="Integer indicating the physics process. Default is 0.")
    parser.add_argument("-g", "--geometry", type=str, default="Ellipse",  help="String indicating the geometry of the PAIReD jet. Default is Ellipse.")
    args = parser.parse_args()

    processDataset(args.dataset, args.outputpath,
                   storageprefix=args.storageprefix,
                   filename=args.filename,
                   batchsize=args.batchsize,
                   first_file=args.start_file,
                   N_files=args.N_files,
                   N_update=args.N_update,
                   N_events_to_process=args.N_events,
                   N_threads=args.N_threads,
                   physics_process=args.physicsprocess,
                   PAIReD_geometry=args.geometry)
