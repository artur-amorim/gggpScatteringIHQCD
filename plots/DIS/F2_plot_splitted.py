import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load MC predicted data
F2_pred = pd.read_csv("F2_pred_2.txt", sep = '\t')
F2_pred_Q2max_10 = F2_pred[F2_pred.Q2 < 8.5]
F2_pred_Q2min_10 = F2_pred[F2_pred.Q2 >= 8.5]

# Load experimental data
F2_exp = pd.read_csv("F2_data.txt", sep = '\t')

# Start the plot
fig, ax = plt.subplots(1,2,sharex=False, sharey=False, figsize = (12,5))
x_min = 10**(-6)
x_max = 0.01
ax[0].set_xlim([x_min, x_max])
ax[0].set_xlabel("x")
ax[0].set_xscale("log")
colours = ["r","b","g"]
Q2s = F2_pred_Q2max_10.Q2.unique()
n = len(Q2s)
for i in range(n) :
    # plot the predicted points
    x_values = np.array(F2_pred_Q2max_10[F2_pred_Q2max_10.Q2 == Q2s[i]].x)
    struct_func_values = np.array(F2_pred_Q2max_10[F2_pred_Q2max_10.Q2 == Q2s[i]].F2)
    ax[0].plot(x_values, struct_func_values, colours[i%3], linewidth = 0.5)
    # plot the experimental points
    x_values = np.array(F2_exp[F2_exp.Q2 == Q2s[i]].x)
    struct_func_values = np.array(F2_exp[F2_exp.Q2 == Q2s[i]].F2)
    struct_func_errors = np.array(F2_exp[F2_exp.Q2 == Q2s[i]].F2err)
    pt_fmt = colours[i%3] + "o"
    ax[0].errorbar(x_values, struct_func_values, struct_func_errors, fmt=pt_fmt, markersize=1)
    x_annotate = 0.5 * x_values.min()
    y_annotate = 1.1 * struct_func_values.max()
    text = ax[0].annotate(str(Q2s[i]), xy = (x_annotate, y_annotate),
                                                    xytext=(x_annotate, y_annotate), color = colours[i%3])
    text.set_fontsize(6)
ax[0].set_ylim([0.1, 1.4])

x_min = 1.5*10**(-4)
ax[1].set_xlim([x_min, x_max])
ax[1].set_xlabel("x")
ax[1].set_xscale("log")
colours = ["r","b","g"]
Q2s = F2_pred_Q2min_10.Q2.unique()
n = len(Q2s)
for i in range(n) :
    # plot the predicted points
    x_values = np.array(F2_pred_Q2min_10[(F2_pred_Q2min_10.Q2 == Q2s[i]) & (F2_pred_Q2min_10.x >= x_min)].x)
    struct_func_values = np.array(F2_pred_Q2min_10[(F2_pred_Q2min_10.Q2 == Q2s[i]) & (F2_pred_Q2min_10.x >= x_min)].F2)
    ax[1].plot(x_values, struct_func_values, colours[i%3], linewidth = 0.5)
    # plot the experimental points
    x_values = np.array(F2_exp[F2_exp.Q2 == Q2s[i]].x)
    struct_func_values = np.array(F2_exp[F2_exp.Q2 == Q2s[i]].F2)
    struct_func_errors = np.array(F2_exp[F2_exp.Q2 == Q2s[i]].F2err)
    pt_fmt = colours[i%3] + "o"
    ax[1].errorbar(x_values, struct_func_values, struct_func_errors, fmt=pt_fmt, markersize=1)
    x_annotate = 2.0 * 10**(-3)
    y_annotate = np.min(F2_pred_Q2min_10[(F2_pred_Q2min_10.Q2 == Q2s[i]) & (F2_pred_Q2min_10.x - x_annotate < 1.5*10**(-4))].F2)
    text = ax[1].annotate(str(Q2s[i]), xy = (x_annotate, y_annotate),
                                                    xytext=(x_annotate, y_annotate), color = colours[i%3])
    text.set_fontsize(6)
ax[1].set_ylim(0.6, 1.75)

fig.tight_layout()
fig.suptitle(r'$F_2(x,Q^2)$ vs x')
fig.subplots_adjust(top=0.9)
fig.savefig("Holographic_F2_splitted_2.pdf", dpi = 100)
