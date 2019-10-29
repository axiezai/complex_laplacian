"""
Use the basinhopping algorithm to find best alpha, speed, and frequency
that produces the best DICE scores for a given canonical network
"""

# number stuff imports
import h5py
import numpy as np
import pandas as pd
from scipy.optimize import basinhopping
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

import sys
import os
import time

# spectrome imports
from spectrome.brain import Brain
from spectrome.utils import functions, path
from spectrome.forward import eigenmode

# Limit number of threads
# os.environ["OMP_NUM_THREADS"] = "2"
# os.environ["MKL_NUM_THREADS"] = "2"
# os.environ["NUMEXPR_NUM_THREADS"] = "2"

# hcp template connectome directory
hcp_dir = "../data"

HCP_brain = Brain.Brain()
HCP_brain.add_connectome(hcp_dir)
HCP_brain.reorder_connectome(HCP_brain.connectome, HCP_brain.distance_matrix)
HCP_brain.bi_symmetric_c()
HCP_brain.reduce_extreme_dir()

# Load Pablo's Yeo 2017 canonical network maps
com_dk = np.load("../data/com_dk.npy", allow_pickle=True).item()
DK_df_normalized = pd.read_csv("../data/DK_dictionary_normalized.csv").set_index(
    "Unnamed: 0"
)

# binarize:
ub, lb = 1, 0  # define binary boundaries

DKfc_binarized = pd.DataFrame(
    [], index=DK_df_normalized.index, columns=DK_df_normalized.columns
)
for name in DK_df_normalized.index:
    u = np.mean(np.nan_to_num(DK_df_normalized.loc[name].values))
    s = np.std(np.nan_to_num(DK_df_normalized.loc[name].values))
    threshold = u - s * 0.1
    DKfc_binarized.loc[name] = np.where(
        DK_df_normalized.loc[name].values > threshold, ub, lb
    )


def laplacian_dice(x, Brain, FC_networks, network_name):
    #    start = time.time()
    # Define frequency of interest
    w = 2 * np.pi * x[0]

    # Laplacian, Brain already prep-ed with connectomes outside of function:
    Brain.add_laplacian_eigenmodes(w=w, alpha=x[1], speed=x[2], num_ev=86)

    # binarize per canonical network:
    thresh_vec = np.linspace(0.1, 0.8, 30)
    binary_count = np.zeros(thresh_vec.shape)
    for i in np.arange(0, len(thresh_vec)):
        binary_mat = np.where(Brain.norm_eigenmodes > thresh_vec[i], 1, 0)
        binary_count[i] = np.count_nonzero(binary_mat)

    num_canonical = np.count_nonzero(FC_networks.loc[network_name].values)
    bin_num = np.abs(binary_count - num_canonical).argmin()
    Brain.binary_eigenmodes = np.where(
        Brain.norm_eigenmodes > thresh_vec[bin_num], 1, 0
    )
    # Brain.binary_eigenmodes = np.where(Brain.norm_eigenmodes > 0.6, 1, 0)

    # Dice
    hcp_dice = eigenmode.get_dice_df(Brain.binary_eigenmodes, FC_networks)
    # Compute mean Dice for chosen network:
    ntw_dice = np.round(hcp_dice[network_name].values.astype(np.double), 3)
    min_dice = np.min(ntw_dice)
    #    end = time.time()
    #    print(end - start)
    return min_dice


def laplacian_corr(x, Brain, FC_networks, network_name):
    # start = time.time()
    # w = 2 * np.pi * x[0]

    # Laplacian, Brain already prep-ed with connectomes outside of function:
    Brain.decompose_complex_laplacian(alpha=x[0], k=x[1], num_ev=86)
    canon_network = np.nan_to_num(FC_networks.loc[network_name].values)

    # compute max correlation for optimization
    corrs = np.zeros([Brain.norm_eigenmodes.shape[1], 1])
    for e in np.arange(0, len(corrs)):
        corrs[e] = -spearmanr(np.squeeze(canon_network), Brain.norm_eigenmodes[:, e])[0]

    max_corr = np.min(corrs)
    # end = time.time()
    # print(end - start)
    return max_corr


class BH_bounds(object):
    def __init__(self, xmax=[5, 600], xmin=[0, 0.1]):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)

    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin


allx0 = np.array(
    [
        [0.5, 5],
        [1, 100],
        [0.8, 50],
        [0.8, 200],
        [0.5, 400],
        [3, 15],
        [5, 250],
        [2, 150],
        [2, 300],
        [1, 500],
    ]
)

bnds = BH_bounds()
print(
    "Starting optimization for {} initial condition {}".format(
        str(sys.argv[1]), str(sys.argv[2])
    )
)

if str(sys.argv[3]) == "dice":
    opt_res = basinhopping(
        laplacian_dice,
        x0=allx0[int(sys.argv[2]), :],
        minimizer_kwargs={"args": (HCP_brain, DKfc_binarized, str(sys.argv[1]))},
        niter=1500,
        T=0.1,
        stepsize=2,
        accept_test=bnds,
        seed=24,
        niter_success=100,
        disp=False,
    )
elif str(sys.argv[3]) == "corr":
    opt_res = basinhopping(
        laplacian_corr,
        x0=allx0[int(sys.argv[2]), :],
        minimizer_kwargs={"args": (HCP_brain, DK_df_normalized, str(sys.argv[1]))},
        niter=1500,
        T=0.1,
        stepsize=2,
        accept_test=bnds,
        seed=24,
        niter_success=100,
        disp=True,
    )

opt_alpha = opt_res["x"][0]
opt_phi = opt_res["x"][1]

# print('optimized output: {}'.format(opt_res))
# Recreate the forward solution:
# w_opt = 2 * np.pi * opt_freq
HCP_brain.decompose_complex_laplacian(alpha=opt_alpha, k=opt_phi)

if str(sys.argv[3]) == "dice":
    # binarize per canonical network:
    thresh_vec = np.linspace(0.1, 0.8, 30)
    binary_count = np.zeros(thresh_vec.shape)
    for i in np.arange(0, len(thresh_vec)):
        binary_mat = np.where(HCP_brain.norm_eigenmodes > thresh_vec[i], 1, 0)
        binary_count[i] = np.count_nonzero(binary_mat)

    num_canonical = np.count_nonzero(DKfc_binarized.loc[str(sys.argv[1])].values)
    bin_num = np.abs(binary_count - num_canonical).argmin()
    HCP_brain.binary_eigenmodes = np.where(
        HCP_brain.norm_eigenmodes > thresh_vec[bin_num], 1, 0
    )
    opt_dice = eigenmode.get_dice_df(HCP_brain.binary_eigenmodes, DKfc_binarized)
    ntw_opt_dice = np.round(opt_dice[str(sys.argv[1])].values.astype(np.double), 3)
    # print('all dice scores: {}'.format(ntw_opt_dice))
    min_opt_dice = np.min(ntw_opt_dice)
    ordered_dice = np.argsort(ntw_opt_dice)
    assert ntw_opt_dice[ordered_dice[1]] >= ntw_opt_dice[ordered_dice[0]]
    # print('Basin hopping final dice: {}'.format(opt_res['fun']))
    # print('Forward calculated dice : {}'.format(min_opt_dice))
    assert min_opt_dice == np.round(opt_res["fun"], 3)
elif str(sys.argv[3]) == "corr":
    canon_network = np.nan_to_num(DK_df_normalized.loc[str(sys.argv[1])].values)
    # compute max correlation for optimization
    corrs = np.squeeze(np.zeros([HCP_brain.norm_eigenmodes.shape[1], 1]))
    for e in np.arange(0, len(corrs)):
        prcorr = spearmanr(np.squeeze(canon_network), HCP_brain.norm_eigenmodes[:, e])[
            0
        ]
        corrs[e] = prcorr
        # print(prcorr)

    ntw_opt_corr = np.round(corrs, 3)
    max_opt_corr = np.max(ntw_opt_corr)
    ordered_corr = np.argsort(-ntw_opt_corr)
    # print(ordered_corr)
    print("basinhop:{}".format(opt_res["fun"]))
    print("forward max:{}".format(max_opt_corr))
    assert ntw_opt_corr[ordered_corr[1]] <= ntw_opt_corr[ordered_corr[0]]
    assert max_opt_corr == -np.round(opt_res["fun"], 3)
# Linear Regression for 10 K's and save in a dictionary:
# K = 11
# if str(sys.argv[3]) == 'dice':
#     # create empty list of dicts:
#     LinReg = []
#     keys = ['num','coef','r2score','ordereigs']
#     for k in np.arange(1,K):
#         selected_eigs = HCP_brain.norm_eigenmodes[:,ordered_dice[0:k]]
#         canon_network = np.nan_to_num(DK_df_normalized.loc[str(sys.argv[1])].values).reshape(-1,1)
#         regr = LinearRegression()
#         regr.fit(canon_network, selected_eigs)
#         c = regr.coef_
#         r2 = regr.score(canon_network, selected_eigs)
#         reg_results = {keys[0]:k, keys[1]:c, keys[2]:r2, keys[3]:ordered_dice[0:k]}
#         LinReg.append(reg_results)
#         print('For K = {}, chosen eigs: {}, coefficients: {} , residual error: {}'.format(k, ordered_dice[0:k], c, r2))

#     opt_res['LinRegResults'] = LinReg

#     file_name = str(sys.argv[1]) + str(sys.argv[2]) + "_BH_dice.h5"
#     file_path = os.path.join(hcp_dir, file_name)
#     path.save_hdf5(file_path, opt_res)
#     print("Optimal result: " , opt_res['x'])
# elif str(sys.argv[3]) == 'corr':
#     # create empty list of dicts:
#     LinReg = []
#     keys = ['num','coef','r2score','ordereigs']
#     for k in np.arange(1,K):
#         selected_eigs = HCP_brain.norm_eigenmodes[:,ordered_corr[0:k]]
#         canon_network = np.nan_to_num(DK_df_normalized.loc[str(sys.argv[1])].values).reshape(-1,1)
#         regr = LinearRegression()
#         regr.fit(canon_network, selected_eigs)
#         c = regr.coef_
#         r2 = regr.score(canon_network, selected_eigs)
#         reg_results = {keys[0]:k, keys[1]:c, keys[2]:r2, keys[3]:ordered_corr[0:k]}
#         LinReg.append(reg_results)
#         print('For K = {}, chosen eigs: {}, coefficients: {} , residual error: {}'.format(k, ordered_corr[0:k], c, r2))

#     opt_res['LinRegResults'] = LinReg

file_name = str(sys.argv[1]) + str(sys.argv[2]) + "_BH_corr.h5"
file_path = os.path.join(hcp_dir, file_name)
path.save_hdf5(file_path, opt_res)
print("Optimal result: ", opt_res["x"])
