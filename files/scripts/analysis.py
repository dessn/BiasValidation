import numpy as np
from chainconsumer import ChainConsumer
from matplotlib import pyplot as plt
from pathlib import Path
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from sklearn import gaussian_process
from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel
import scipy.integrate as integrate


def get_wfit_dirs(pippin_output):
    wfit_dir = Path(pippin_output) / "8_COSMOFIT" / "WFIT"
    wfit_dirs = list(wfit_dir.iterdir())
    return wfit_dirs

def get_wfit_files(wfit_dir):
    out_dir = wfit_dir / "output" / "CHI2GRID"
    wfit_files = list(out_dir.iterdir())
    return wfit_files

def read_wfit(wfit_file):
    chi2 = np.loadtxt(wfit_file, dtype=str)
    Om = np.array([float(i[3]) for i in chi2[1:]])
    w0 = np.array([float(i[2]) for i in chi2[1:]])
    weight = np.array([float(i[5]) for i in chi2[1:]])
    weight = np.exp(-0.5 * weight)
    ind = np.argmax(weight)
    best_Om = Om[ind]
    best_w0 = w0[ind]
    best_weight = weight[ind]
    data = [Om, w0, weight]
    best = [best_Om, best_w0, best_weight]
    return wfit_file, data, best

def unpack_dict(wfit_dict):
    Om_list = np.hstack([l[0] for l in [wfit_dict[k]["data"] for k in wfit_dict.keys()]])
    w0_list = np.hstack([l[1] for l in [wfit_dict[k]["data"] for k in wfit_dict.keys()]])
    weight_list = np.hstack([l[2] for l in [wfit_dict[k]["data"] for k in wfit_dict.keys()]])
    best_Om_list = np.hstack([l[0] for l in [wfit_dict[k]["best"] for k in wfit_dict.keys()]])
    best_w0_list = np.hstack([l[1] for l in [wfit_dict[k]["best"] for k in wfit_dict.keys()]])
    best_weight_list = np.hstack([l[2] for l in [wfit_dict[k]["best"] for k in wfit_dict.keys()]])
    return Om_list, w0_list, weight_list, best_Om_list, best_w0_list, best_weight_list

def get_closest(wfit_dict):
    Om_list, w0_list, weight_list, best_Om_list, best_w0_list, best_weight_list = unpack_dict(wfit_dict)
    Om = np.mean(best_Om_list)
    Om_std = np.std(best_Om_list)
    w0 = np.mean(best_w0_list)
    w0_std = np.std(best_w0_list)
    dist = np.sqrt(np.power(best_Om_list - Om, 2) + np.power(best_w0_list - w0, 2))
    ind = np.argmin(dist)
    best_file = list(wfit_dict.keys())[ind]
    return read_wfit(best_file)

def plot_nominal(wfits, options, output):
    print("Plotting nominal")
    f, data, best = get_closest(wfits["FR"][list(wfits["FR"].keys())[0]])
    c = ChainConsumer()
    c.add_chain([data[0], data[1]], name="Nominal Cosmology", parameters=["Om", "w0"], weights=data[2], grid=True, shade_alpha=0.2)
    fig = c.plotter.plot(figsize=(12, 8))
    ax = fig.get_axes()[2]
    ax.scatter([best[0]], [best[1]], label="Nominal Cosmology Input", s=1)
    for (i, k) in enumerate(wfits["V"].keys()):
        key = list(wfits["V"][k].keys())[0]
        Om = wfits["V"][k][key]["Om"]
        w0 = wfits["V"][k][key]["w0"]
        ax.scatter([Om], [w0], label=f"Validation Cosmology {i+1} Input", s=1)
    Om = options.get("Om", [])
    w0 = options.get("w0", [])
    for i in range(len(Om)):
        ax.scatter([Om[i]], [w0[i]], label=f"Proposal Cosmology {i+1} Input", s=1)
    ax.legend()
    plt.savefig(output / "nominal.svg")
    plt.close()

def get_GP(wfit_dict):
    Om_list, w0_list, weight_list, best_Om_list, best_w0_list, best_weight_list = unpack_dict(wfit_dict)
    X = best_Om_list.reshape(-1, 1)
    Y = best_w0_list
    kernel = ConstantKernel() + Matern(length_scale=2, nu=3/2) + WhiteKernel(noise_level=1)
    gp = gaussian_process.GaussianProcessRegressor(kernel=kernel)
    gp.fit(X, Y)
    return gp

def apply_GP(wfits):
    rtn = {"FR": {}, "V": {}}
    k = list(wfits["FR"].keys())[0]
    wfit_dict = wfits["FR"][k]
    Om_list, w0_list, weight_list, best_Om_list, best_w0_list, best_weight_list = unpack_dict(wfit_dict)
    mean_Om = np.mean(best_Om_list)
    gp = get_GP(wfit_dict)
    w0_correction = gp.predict(best_Om_list.reshape(-1,1), return_std=False)
    best_w0_list -= w0_correction
    best_Om_list -= mean_Om
    rtn["FR"][k] = (best_Om_list, best_w0_list)
    for (i, k) in enumerate(wfits["V"].keys()):
        wfit_dict = wfits["V"][k]
        Om_list, w0_list, weight_list, best_Om_list, best_w0_list, best_weight_list = unpack_dict(wfit_dict)
        mean_Om = np.mean(best_Om_list)
        gp = get_GP(wfit_dict)
        w0_correction = gp.predict(best_Om_list.reshape(-1,1), return_std=False)
        best_w0_list -= w0_correction
        best_Om_list -= mean_Om
        rtn["V"][k] = (best_Om_list, best_w0_list)
    return rtn

def plot_GPE(wfits, options, output):
    print("Plotting GPE")
    k = list(wfits["FR"].keys())[0]
    wfit_dict = wfits["FR"][k]
    gp = get_GP(wfit_dict)
    Om_list, w0_list, weight_list, best_Om_list, best_w0_list, best_weight_list = unpack_dict(wfit_dict)
    x = np.linspace(best_Om_list.min(), best_Om_list.max(), 100)
    y, dy = gp.predict(x.reshape(-1, 1), return_std=True)
    plt.scatter(best_Om_list, best_w0_list, label="Nominal Cosmology", s=1)
    plt.errorbar(x, y, yerr=dy)
    for (i, k) in enumerate(wfits["V"].keys()):
        wfit_dict = wfits["V"][k]
        gp = get_GP(wfit_dict)
        Om_list, w0_list, weight_list, best_Om_list, best_w0_list, best_weight_list = unpack_dict(wfit_dict)
        x = np.linspace(best_Om_list.min(), best_Om_list.max(), 100)
        y, dy = gp.predict(x.reshape(-1, 1), return_std=True)
        plt.scatter(best_Om_list, best_w0_list, label=f"Validation Cosmology {i+1}", s=1)
        plt.errorbar(x, y, yerr=dy)
    plt.title("Gaussian Process Fit")
    plt.xlabel("Om")
    plt.ylabel("w0")
    plt.legend()
    plt.savefig(output / "GPE.svg")
    plt.close()

def min_kde(Om, w0, bandwidth):
    scale = 10.0
    Om_min = min(Om)
    Om_min -= scale * abs(Om_min)
    Om_max = max(Om) 
    Om_max += scale * abs(Om_max)
    w0_min = min(w0) 
    w0_min -= scale * abs(w0_min)
    w0_max = max(w0) 
    w0_max += scale * abs(w0_max)
    X, Y = np.mgrid[Om_min:Om_max:2000j, w0_min:w0_max:2000j]
    xy_train = np.vstack([w0, Om]).T
    xy_sample = np.vstack([Y.ravel(), X.ravel()]).T

    def _scores(estimator, sample):
        print(sample.shape)
        scores = estimator.score_samples(sample)
        scores = scores[scores != float('-inf')]
        return np.mean(scores)

    kernels = ['exponential', 'gaussian']
    grid = GridSearchCV(KernelDensity(), {'bandwidth': bandwidth, 'kernel': kernels}, scoring=_scores)
    grid.fit(xy_train)
    kde = grid.best_estimator_
    print(f"Kernel: {kde.kernel}, Bandwidth: {kde.bandwidth}")
    dens = np.reshape(np.exp(kde.score_samples(xy_sample)), X.shape)
    return X, Y, dens

def integrate_KDE(X, Y, dens):
    x = X[:, 0]
    y = Y[0, :]
    return integrate.simps(integrate.simps(dens, y), x)

def get_probability(X, Y, dens):
    prob = np.ones(dens.shape)
    x = X.copy()
    y = Y.copy()
    d = dens.copy()
    zs = np.linspace(dens.min(), dens.max(), 200)
    for z in zs:
        mask = dens < z
        d[mask] = 0
        p = integrate_KDE(x, y, d)
        prob[~mask] = p
    return prob

def calculate_KDE(wfits):
    rtn = {"FR": {}, "V": {}}
    gp = apply_GP(wfits)
    k = list(wfits["FR"].keys())[0]
    Om, w0 = gp["FR"][k]
    bandwidth = np.linspace(0.02, 0.05, 100)
    X, Y, Z = min_kde(Om, w0, bandwidth)
    rtn["FR"][k] = (X, Y, Z)
    for (i, k) in enumerate(wfits["V"].keys()):
        Om, w0 = gp["V"][k]
        X, Y, Z = min_kde(Om, w0, bandwidth)
        rtn["V"][k] = (X, Y, Z)
    return rtn

def plot_KDE(wfits, options, output):
    print("Plotting KDE")
    percentile = [0.6827, 0.9545]
    #percentile = [1 - 0.6827, 1 - 0.9545]
    #percentile.reverse()
    def fmt(x):
        x = 100 * x 
        s = f"{x:.4f}"
        if s.endswith("0"):
            s = f"{x:.0f}"
        return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"
    
    gp = apply_GP(wfits)
    kde = calculate_KDE(wfits)

    n = int(np.ceil(np.sqrt(len(gp["V"]))))
    fig, axes = plt.subplots(n, n, layout="constrained")
    for (i, k) in enumerate(gp["V"].keys()):
        x = int(i % n)
        y = int(np.floor(i / n))
        ax = axes[x, y]
        Om, w0 = gp["V"][k]
        X, Y, Z = kde["V"][k]
        print(f"Integral: {integrate_KDE(X, Y, Z)}")
        Z = get_probability(X, Y, Z)
        ax.scatter(Om, w0, s=1)
        ax.set_title(f"Validation Cosmology {i + 1}")
        CS = ax.contour(X, Y, Z, levels=percentile, linewidths=1, linestyles='dashed') 
        plt.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=8)
        ax.set_xlabel("Om")
        ax.set_ylabel("w0 - GPE")
        scale = 1.0
        Om_min = min(Om)
        Om_min -= scale * abs(Om_min)
        Om_max = max(Om) 
        Om_max += scale * abs(Om_max)
        w0_min = min(w0) 
        w0_min -= scale * abs(w0_min)
        w0_max = max(w0) 
        w0_max += scale * abs(w0_max)
        ax.set_xlim(Om_min, Om_max)
        ax.set_ylim(w0_min, w0_max)

    fig.suptitle("KDE Contours")
    plt.savefig(output / "KDE.svg")
    plt.close()
