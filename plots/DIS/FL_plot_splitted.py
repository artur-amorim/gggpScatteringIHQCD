import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load MC predicted data
FL_pred = pd.read_csv("FL_pred.txt", sep = '\t')
FL_pred_Q2max_10 = FL_pred[FL_pred.Q2 <= 10]
FL_pred_Q2min_10 = FL_pred[FL_pred.Q2 > 10]

# Load experimental data
FL_exp = pd.read_csv("FL_data.txt", sep = '\t')

# Find unique values of Q^2 in FL data with Q2 <= 10
Q2s_max_10 = FL_pred_Q2max_10.Q2.unique()
n_Q2max_10 = len(Q2s_max_10)

# Find unique values of Q^2 in FL data with Q2 > 10
Q2s_min_10 = FL_pred_Q2min_10.Q2.unique()
n_Q2min_10 = len(Q2s_min_10)

# number of rows
n_rows = max(n_Q2max_10, n_Q2min_10)

# Start the plot
fig, ax = plt.subplots(n_rows,2,sharex=False, sharey=False, figsize = (12,10))
fig.suptitle(r'$F_L(x,Q^2)$ vs x')

j_Q2max_10, j_Q2min_10 = 0, 0

for i in range(n_rows * 2) :
    if i == (n_rows*2-1):
        break
    figp = plt.subplot(n_rows,2,i+1)
    if (i % 2 == 1 and j_Q2min_10 < n_Q2min_10):
        x_min = 2*10**(-4)
        x_max = 3*10**(-3)
        figp.set_xlabel("x")
        x_values = np.array(FL_pred_Q2min_10[FL_pred_Q2min_10.Q2 == Q2s_min_10[j_Q2min_10]].x)
        struct_func_values = np.array(FL_pred_Q2min_10[FL_pred_Q2min_10.Q2 == Q2s_min_10[j_Q2min_10]].FL)
        plt.plot(x_values, struct_func_values, "b", linewidth = 0.5)
        plt.plot(x_values,0*x_values,"k--")
        x_annotate = 1.5*10**(-3)
        y_annotate = 0.5*(struct_func_values.min()+struct_func_values.max())
        text = plt.annotate(r'$Q^2 = ' + str(Q2s_min_10[j_Q2min_10]) + '$', xy = (x_annotate, y_annotate),
                                                    xytext=(1.1 * x_annotate, 0.7 * y_annotate))
        text.set_fontsize(10)
        # plot the experimental points
        x_values = np.array(FL_exp[FL_exp.Q2 == Q2s_min_10[j_Q2min_10]].x)
        struct_func_values = np.array(FL_exp[FL_exp.Q2 == Q2s_min_10[j_Q2min_10]].FL)
        struct_func_errors = np.array(FL_exp[FL_exp.Q2 == Q2s_min_10[j_Q2min_10]].FLerr)
        plt.errorbar(x_values, struct_func_values, struct_func_errors, fmt="bo", markersize=1)
        plt.xscale("log")
        plt.xlim(x_min, x_max)
        j_Q2min_10 += 1
    elif (i % 2 == 0 and j_Q2max_10 < n_Q2max_10):
        x_min = 2*10**(-5)
        x_max = 8*10**(-4)
        figp.set_xlabel("x")
        x_values = np.array(FL_pred_Q2max_10[FL_pred_Q2max_10.Q2 == Q2s_max_10[j_Q2max_10]].x)
        struct_func_values = np.array(FL_pred_Q2max_10[FL_pred_Q2max_10.Q2 == Q2s_max_10[j_Q2max_10]].FL)
        plt.plot(x_values, struct_func_values, "b",linewidth = 0.5)
        plt.plot(x_values,0*x_values,"k--")
        x_annotate = 4*10**(-4)
        y_annotate = 0.5*(struct_func_values.min()+struct_func_values.max())
        text = plt.annotate(r'$Q^2 = ' + str(Q2s_max_10[j_Q2max_10]) + '$', xy = (x_annotate, y_annotate),
                                                    xytext=(1.1 * x_annotate, 0.7 * y_annotate))
        text.set_fontsize(10)
        # plot the experimental points
        x_values = np.array(FL_exp[FL_exp.Q2 == Q2s_max_10[j_Q2max_10]].x)
        struct_func_values = np.array(FL_exp[FL_exp.Q2 == Q2s_max_10[j_Q2max_10]].FL)
        struct_func_errors = np.array(FL_exp[FL_exp.Q2 == Q2s_max_10[j_Q2max_10]].FLerr)
        plt.errorbar(x_values, struct_func_values, struct_func_errors, fmt="bo", markersize=1)
        plt.xscale("log")
        plt.xlim(x_min, x_max)
        j_Q2max_10 += 1

fig.delaxes(ax[n_rows-1][1])
l  = 0.05  # the left side of the subplots of the figure
r = 0.975    # the right side of the subplots of the figure
b = 0.05   # the bottom of the subplots of the figure
t = 0.95      # the top of the subplots of the figure
ws = 0.1   # the amount of width reserved for blank space between subplots
hs = 0.3  # the amount of height reserved for white space between subplots
plt.subplots_adjust(left = l, right = r, top = t, bottom = b, wspace = ws, hspace = hs)
#fig.tight_layout()
plt.savefig("Holographic_FL_splitted.pdf", dpi = 100)
plt.close()
