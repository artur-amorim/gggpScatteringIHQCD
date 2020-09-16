# Import pandas, matplotlib and numpy
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

# import sys to read CMD input
import sys

def plotPDFsNLO(Q2list, xg_ihqcd_path, coupling = ""):
    """
        Given a list of Q2 values, it plots and saves the
        plots that compare the PDFs given by CT18 and NNPDF NLO
        with the ones predicted by IHQCD
    """
    gluonPDFIHQCD = pd.read_csv(xg_ihqcd_path, sep = '\t')
    for Q2 in Q2list:
        # Get data and the errors from NNPDF
        gluonPDFDataNNPDF = pd.read_csv("NNPDF31_nlo_as_0118_nf_4_Q2_"+str(Q2)+".txt", sep = '\t')
        # Get data and the errors from CT18NNLO
        gluonPDFDataCT18 = pd.read_csv("CT18NLO_Q2_"+str(Q2)+".txt", sep = '\t')
        # Get IHQCD prediction for the gluon pdfs
        gluonPDFIHQCDQ2 = gluonPDFIHQCD[gluonPDFIHQCD.Q2 == Q2]
        # Computing h(Q^2)
        holographic_xg = float(gluonPDFIHQCDQ2[gluonPDFIHQCDQ2.x == 0.01].xg)
        xg_CT18 = float(gluonPDFDataCT18[gluonPDFDataCT18.x==0.01].xg)
        xg_NNPDF = float(gluonPDFDataNNPDF[gluonPDFDataNNPDF.x==0.01].xg)
        avg_qcd_xg = 0.5 * (xg_CT18 + xg_NNPDF)
        h_Q2 = (avg_qcd_xg - holographic_xg)
        xs = np.array(gluonPDFIHQCDQ2.x)
        xg = np.array(gluonPDFIHQCDQ2.xg)
        xg_error = np.array(gluonPDFIHQCDQ2.xgError)
        xg = xg + h_Q2
        plt.figure()
        # Plot the data from NNPDF
        plt.plot(gluonPDFDataNNPDF.x, gluonPDFDataNNPDF.xg, "b-", label = r'NNPDF $Q^2 = $'+str(Q2))
        # Plot the data from CT18NNLO
        plt.plot(gluonPDFDataCT18.x, gluonPDFDataCT18.xg, "r-", label = r'CT18 $Q^2 = $'+str(Q2))
        # Plt the data from IHQCD
        plt.plot(xs, xg, "k-", linewidth = 0.1, label = r'IHQCD $Q^2 = $'+str(Q2))
        plt.xlabel("x")
        plt.xscale("log")
        plt.xlim(10**(-6),0.01)
        plt.ylabel(r'$xg(x,Q^2)$')
        plt.fill_between(gluonPDFDataNNPDF.x, gluonPDFDataNNPDF.xg - gluonPDFDataNNPDF["error-"],
                 gluonPDFDataNNPDF.xg + gluonPDFDataNNPDF["error+"], facecolor='blue', alpha = 0.2)
        plt.fill_between(gluonPDFDataCT18.x, gluonPDFDataCT18.xg - gluonPDFDataCT18["error-"],
                 gluonPDFDataCT18.xg + gluonPDFDataCT18["error+"], facecolor='red', alpha = 0.2)
        plt.fill_between(xs, xg - xg_error, xg + xg_error, facecolor = 'yellow', alpha = 1)
        plt.legend()
        plt.savefig("gluon_PDF_Q2_"+str(Q2)+ "_"+coupling + "_NLO.pdf")
        plt.show()


# Read input from command line
argument_vector = sys.argv
if(len(argument_vector) == 1):
    print("Usage: " + argument_vector[0] + " Q2 vals.")
else:
    # The first value of argument_vector is the program name so we exclude it
    Q2vals = list(map(int,argument_vector[1:]))
    # Do the plots now
    #plotPDFsNNLO(Q2vals, "with_uncertainty")
    plotPDFsNLO(Q2vals, "xg_IHQCD_pred.txt")
