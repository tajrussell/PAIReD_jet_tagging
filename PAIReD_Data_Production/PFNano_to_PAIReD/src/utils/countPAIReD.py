"""count PAIReD jets

This script allows the user to count the PAIReD jets in all files of a
directory. It saves the information in a text file in the respective
directory:
    PAIReD__counts.txt

This tool accepts ROOT files (.root) in the format of PAIReD files.

This script requires the following modules to be installed within the
Python environment you are running this script in:

    * uproot
    * awkward
    * glob
"""



# import all modules needed
import argparse
import uproot
import awkward as ak
import numpy as np
import glob



def group(n, sep = ','):
    s = str(abs(int(n)))[::-1]
    groups = []
    i = 0
    while i < len(s):
        groups.append(s[i:i+3])
        i+=3
    retval = sep.join(groups)[::-1]
    if n < 0:
        return '-%s' % retval
    else:
        return retval

def group_list(list, sep = ','):
    out = "("
    for n in list[:-1]:
        out += "%10s" % group(n, sep)
        out += " | "
    out += "%10s" % group(list[-1], sep)
    out += ")"
    return out

# ******************************************************************************


"""
 * Function: countPAIReD  [main function]
 * --------------------------------------
Count PAIReD jets in directory.

Parameters
----------
directoryPath : str
    Paths to the directory (PAIReD training file)

Returns
-------

    nothing, but creates output file with counts
"""

def countPAIReD(directoryPath, labels=["CC", "BB", "ll", "cl", "cc", "bl", "bb"], flavs=[11,13,15,0]):

    print("\n*****************************************************************")
    print("  Start countPAIReD()")
    print("*****************************************************************\n")

    print("** Look at directory:", directoryPath)
    print("** Consider the labels:", labels)
    print("** and vector boson flavors:", flavs)

    txt = "Counts of PAIReD jets in the training set:\n" + directoryPath
    txt += "\n\ndecomposed in labels:\n" + str(labels)
    txt += "\nand vector boson flavors:\n" + str(flavs)

    # save current cursor position for later
    pos = len(txt)

    txt += "\n\n################################################"
    txt += "\n   Numbers of individual files"
    txt += "\n################################################"

    # get all root files in directory
    files = glob.glob(directoryPath + "/*PAIReD*.root")

    # total number of PAIReD jets accross all files
    Ntot = {l : np.zeros(len(flavs)) for l in labels}

    # total number of PAIReD jets for testing (test set)
    Ntest = {l : np.zeros(len(flavs)) for l in labels}

    # total number of inactive PAIReD jets (inactive)
    Ninactive = {l : np.zeros(len(flavs)) for l in labels}

    # total number of training PAIReD jets (training set)
    Ntrain = {l : np.zeros(len(flavs)) for l in labels}

    print("\n** Loop through all files of the training set")
    # loop trough all files and count labels
    for f in files:
        with uproot.open(f) as inputFile:
            inputTree = inputFile['tree']

            # branches to load
            branches = ["label_%s" % l for l in labels] + ["MC_vector_flav"]

            # load labels
            Events = inputTree.arrays(branches)

            N = {l : np.zeros(len(flavs)) for l in labels}
            Nsum = np.zeros(len(flavs))

            # number of PAIReD jets in this file
            for l in labels:
                # count flavs for given label in file
                Nl = []
                for flav in flavs:
                    Nl.append(ak.sum((Events["label_%s"%l]==1) & (Events.MC_vector_flav==flav)))
                N[l] = np.array(Nl)
                Nsum += N[l]

                # update Ntot, Ntest, Ninactive, Ntrain
                Ntot[l] += N[l]
                if "test_PAIReD" in f:
                    Ntest[l] += N[l]
                elif "inactive_PAIReD" in f:
                    Ninactive[l] += N[l]
                else:
                    Ntrain[l] += N[l]

            # Write to txt file
            txt += "\n\nfile: " + f
            txt += ("\n** Total number of PAIReD jets in file: \n   %s" % group(np.sum(Nsum)) +
                    "\n   decomposed in vector flavs %s: \n     %s" % (str(flavs), group_list(Nsum)))
            for l in labels:
                txt += ("\n** %s PAIReD jets:\n   %s" % (l.upper(), group(np.sum(N[l]))) +
                        "\n     %s" % (group_list(N[l])))
    
    # total of each category
    for N in [Ntot, Ntest, Ninactive, Ntrain]:
        N["tot"] = np.zeros(len(flavs))
        for l in labels:
            N["tot"] += N[l]
    
    # calculate test set percentages
    p_test = {}
    for l in labels+["tot"]:
        p_test[l] = np.sum(Ntest[l])/np.sum(Ntot[l]-Ninactive[l])*100 if np.sum(Ntot[l]-Ninactive[l])>0 else 0
    
    # calculate inactive percentages
    p_inactive = {}
    for l in labels+["tot"]:
        p_inactive[l] = np.sum(Ninactive[l])/np.sum(Ntot[l])*100 if np.sum(Ntot[l])>0 else 0
    

    # total numbers
    summarytxt = ("\n\n** Total number of PAIReD jets in data set: \n   %s" % group(np.sum(Ntot["tot"])) +
                  "\n   decomposed in vector flavors %s: \n   %13s  %s" % (str(flavs), "", group_list(Ntot["tot"])))
    for l in labels:
        summarytxt += ("\n** %s PAIReD jets:\n   %13s" % (l.upper(), group(np.sum(Ntot[l]))) +
                "  %s" % (group_list(Ntot[l])))

    # test+train numbers
    summarytxt += ("\n\n** Total number of PAIReD jets in (training+test) set: \n   %s" % group(np.sum(Ntot["tot"] - Ninactive["tot"])) +
                  "\n   decomposed in vector flavors %s: \n   %13s  %s" % (str(flavs), "", group_list(Ntot["tot"] - Ninactive["tot"])))
    for l in labels:
        summarytxt += ("\n** %s PAIReD jets:\n   %13s" % (l.upper(), group(np.sum(Ntot[l] - Ninactive[l]))) +
                "  %s" % (group_list(Ntot[l] - Ninactive[l])))

    # test numbers
    summarytxt += ("\n\n** Total number of PAIReD jets in test set: \n   %s   (%.1f %s)" % (group(np.sum(Ntest["tot"])),p_test["tot"],"%") +
                  "\n   decomposed in vector flavors %s: \n   %13s  %s" % (str(flavs), "", group_list(Ntest["tot"])))
    for l in labels:
        summarytxt += ("\n** %s PAIReD jets:\n   %13s" % (l.upper(), group(np.sum(Ntest[l]))) +
                "  %s" % (group_list(Ntest[l])))
        
    # inactive numbers
    summarytxt += ("\n\n** Total number of PAIReD jets with status *inactive*: \n   %s   (%.1f %s)" % (group(np.sum(Ninactive["tot"])),p_inactive["tot"],"%") +
                  "\n   decomposed in vector flavors %s: \n   %13s  %s" % (str(flavs), "", group_list(Ninactive["tot"])))
    for l in labels:
        summarytxt += ("\n** %s PAIReD jets:\n   %13s" % (l.upper(), group(np.sum(Ninactive[l]))) +
                "  %s" % (group_list(Ninactive[l])))
        
    
    # HTML table
    summarytxt += ('\n\nHTML table cell:\n<td class="tg-ltxa"><b>' +
                  group(np.sum(Ntot["tot"])) + ' (%.0f%s)</b><br>(' %(p_test["tot"],"%"))
    for l in labels[:-1]:
        summarytxt += group(np.sum(Ntot[l])) + ' : '
    
    summarytxt += group(np.sum(Ntot[labels[-1]])) + ')</td>'

    # Excel table
    if flavs==[11,13,15,0] and labels==["ll", "cc", "bb", "cl", "bl", "bb_ttbar"]:
        summarytxt += "\n\nExcel table input:\n"
        for i in range(4):
            for l in labels:
                summarytxt += "%i %i " % (Ntot[l][i] - Ninactive[l][i], Ntest[l][i])
            summarytxt += "\n"


    # create txt output file
    txtfilename = "PAIReD__counts.txt"
    txtfile = open(directoryPath + "/" + txtfilename, "w+")
    txt = txt[:pos] + summarytxt + txt[pos:]
    txtfile.write(txt)
    txtfile.close()
    print("** Created output file for counts:", txtfilename)

    print(summarytxt)

    print("\n*****************************************************************")
    print("  End countPAIReD()")
    print("*****************************************************************\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count PAIReD jets in training set under given directory path.")
    parser.add_argument("dirpath", help="Path to the PAIReD training set (directory)")
    args = parser.parse_args()

    countPAIReD(args.dirpath)
