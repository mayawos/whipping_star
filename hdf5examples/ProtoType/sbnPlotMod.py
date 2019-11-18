import h5py
import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import matplotlib as mpl
plt.style.use('ggplot')


def chi2Fit(data, nbins=100):
    return 1 # remove this to do chi2 calculation
    from scipy.stats import chi2
    from scipy.optimize import curve_fit
    vals, edges = np.histogram(data, bins=nbins)
    center = (edges[:-1] + edges[1:]) / 2
    cfa, cfb = curve_fit(lambda x,y :chi2.pdf(x,y), center, vals)
    return cfa[0]

def chi2Fit2(data):
    return 1 # remove this to do chi2 calculation
    from scipy.stats import chi2
    a,b,c = chi2.fit(data)
    return a



def addContour(x,y,z, val):
    import scipy.interpolate
    N=100

    levels = [1,2,3]

    xi = np.linspace(x.min(), x.max(), N)
    yi = np.linspace(y.min(), y.max(), N)
    zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')

    from matplotlib import pyplot as plt
    cs = plt.contour(xi,yi,zi,levels, colors=["red", "blue", "green"])

if __name__=="__main__":
    import sys
    fin  = sys.argv[1]
    fout = sys.argv[2]
    pct = float(sys.argv[3])
    assert(fin != fout)

    F = h5py.File(fin, mode="r")

    # The entries for specified universe
    X = F["gridx"][:]
    X = X.reshape((X.shape[0],))
    Y = F["gridy"][:]
    Y = Y.reshape((Y.shape[0],))

    B = F["best_grid_point"][:]


    I = F["i_univ"][:]
    D = F["delta_chi"][:]
    NEV = F["n_events"][:]
    # D = F["last_chi_min"][:]

    NPOINT = X.shape[0]
    NUNIV = int(I.shape[0] / NPOINT)
    bb = np.unique(B)
    print(bb)
    # plt.hist(B,bins=NPOINT)
    # plt.show()
    BB = B.reshape((NPOINT, NUNIV))
    u, indices = np.unique(B, return_inverse=True)
    star=u[np.argmax(np.bincount(indices))]
    print(X[star], Y[star])


    GG=D.reshape((NPOINT, NUNIV))
    PCT = np.percentile(GG, pct, axis=1)

    NNNN = NEV.reshape((NPOINT, NUNIV))
    MEANNEV = np.mean(NNNN, axis=1)


    EFFDOF =np.array([chi2Fit(g) for g in GG])
    EFFDOF2 =np.array([chi2Fit2(g) for g in GG])


    dX = np.sort(np.unique(X))
    dY = np.sort(np.unique(Y))
    # PCT.reshape((dX.shape[0], dY.shape[0]))
    # x,y = np.meshgrid(dX, dY)
    # plt.pcolormesh(x, y, PCT)


    from matplotlib import pyplot as plt
    import matplotlib as mpl
    plt.style.use('ggplot')

    f, ax = plt.subplots(1, 4, figsize=(24,6))
    im1 = ax[0].pcolor(PCT.reshape((dX.shape[0], dY.shape[0])), vmin=0, vmax=5)
    bar=f.colorbar(im1, ax=ax[0])
    bar.set_label("$\Delta\chi_\\mathrm{crit}^{%.1f}$"%pct)

    im2 = ax[1].pcolor(EFFDOF.reshape((dX.shape[0], dY.shape[0])), vmin=0, vmax=2.5)
    bar=f.colorbar(im2, ax=ax[1])
    bar.set_label("effective degrees of freedom --- curve fit")

    im3 = ax[2].pcolor(EFFDOF2.reshape((dX.shape[0], dY.shape[0])), vmin=0, vmax=2.5)
    bar=f.colorbar(im3, ax=ax[2])
    bar.set_label("effective degrees of freedom --- dist fit")

    im4 = ax[3].pcolor(np.log10(MEANNEV.reshape((dX.shape[0], dY.shape[0]))))
    bar=f.colorbar(im4, ax=ax[3])
    bar.set_label("log10 mean n evnents")

    ax[0].set_xlabel("$\sin^2$ something")
    ax[0].set_ylabel("$\log_{10}(\Delta m^2)$")

    xlocs = [int(i) for i in ax[1].xaxis.get_majorticklocs()]
    xlab = ["{:.1f}".format(dY[i]) for i in xlocs[:-1]]
    ylocs = [int(i) for i in ax[1].yaxis.get_majorticklocs()]
    ylab = ["{:.1f}".format(dX[i]) for i in ylocs[:-1]]

    ax[0].xaxis.set_ticklabels(xlab)
    ax[1].xaxis.set_ticklabels(xlab)
    ax[2].xaxis.set_ticklabels(xlab)
    ax[3].xaxis.set_ticklabels(xlab)
    ax[0].yaxis.set_ticklabels(ylab)
    ax[1].yaxis.set_ticklabels(ylab)
    ax[2].yaxis.set_ticklabels(ylab)
    ax[3].yaxis.set_ticklabels(ylab)

    f.savefig(fout)


    # CC = getContour(Y, X, Z, float(sys.argv[3]), "contour"+fout)

    # # from IPython import embed
    # # embed()
    # # sys.exit(1)

    # plt.close()
    # plt.clf()
    # plt.style.use('ggplot')
    # fig, ax = plt.subplots(1)
    # ax.set_aspect('equal')
    # ax.scatter(Y,X, c=Z, marker="s")
    # plt.savefig(fout)

    # # GG=D2.reshape((NPOINT, NUNIV))
    # # Z=np.sum(GG, axis=1)

    # from matplotlib import pyplot as plt
    # import matplotlib as mpl
    # plt.close()
    # plt.clf()
    # plt.style.use('ggplot')
    # fig, ax = plt.subplots(1)
    # ax.set_aspect('equal')
    # ax.scatter(Y,X, c=Z, marker="s")
    # plt.savefig("chi_"+fout)

