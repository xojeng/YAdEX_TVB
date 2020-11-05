from tvb.simulator.lab import *
from matplotlib.tri import Triangulation
from pylab import *

cortex = cortex.Cortex(load_default=True)

def multiview(data, suptitle='', figsize=(15, 10), **kwds):

    cs = cortex
    vtx = cs.vertices
    tri = cs.triangles
    rm = cs.region_mapping
    x, y, z = vtx.T
    lh_tri = tri[(rm[tri] < 38).any(axis=1)]
    lh_vtx = vtx[rm < 38]
    lh_x, lh_y, lh_z = lh_vtx.T
    lh_tx, lh_ty, lh_tz = lh_vtx[lh_tri].mean(axis=1).T
    rh_tri = tri[(rm[tri] >= 38).any(axis=1)]
    rh_vtx = vtx[rm < 38]
    rh_x, rh_y, rh_z = rh_vtx.T
    rh_tx, rh_ty, rh_tz = vtx[rh_tri].mean(axis=1).T
    tx, ty, tz = vtx[tri].mean(axis=1).T

    views = {
        'lh-lateral': Triangulation(-x, z, lh_tri[argsort(lh_ty)[::-1]]),
        'lh-medial': Triangulation(x, z, lh_tri[argsort(lh_ty)]),
        'rh-medial': Triangulation(-x, z, rh_tri[argsort(rh_ty)[::-1]]),
        'rh-lateral': Triangulation(x, z, rh_tri[argsort(rh_ty)]),
        'both-superior': Triangulation(y, x, tri[argsort(tz)]),
    }

    def plotview(i, j, k, viewkey, z=None, zmin=None,zmax =None, zthresh=None, suptitle='', shaded=True, cmap=plt.cm.coolwarm, viewlabel=False):
        v = views[viewkey]
        ax = subplot(i, j, k)
        if z is None:
            z = rand(v.x.shape[0])
        if not viewlabel:
            axis('off')
        kwargs = {'shading': 'gouraud'} if shaded else {'edgecolors': 'k', 'linewidth': 0.1}
        if zthresh:
            z = z.copy() * (abs(z) > zthresh)
        tc = ax.tripcolor(v, z, cmap=cmap, **kwargs)
        tc.set_clim(vmin=-zmin, vmax=zmax)
        ax.set_aspect('equal')
        if suptitle:
            ax.set_title(suptitle, fontsize=24)
        if viewlabel:
            xlabel(viewkey)

    figure(figsize=figsize)
    plotview(2, 3, 1, 'lh-lateral', data, **kwds)
    plotview(2, 3, 4, 'lh-medial', data, **kwds)
    plotview(2, 3, 3, 'rh-lateral', data, **kwds)
    plotview(2, 3, 6, 'rh-medial', data, **kwds)
    plotview(1, 3, 2, 'both-superior', data, suptitle=suptitle, **kwds)
    subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0, wspace=0, hspace=0)
    plt.show()

import matplotlib.animation as manimation
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
def animation():
    # example for animation

    begin = 0
    end = 100

    my_dpi=100 # resolution
    # parameter fig
    text_fontsize = 20.0
    title_fontsize = 20.0
    fig = plt.figure(figsize=(2000/my_dpi, 800/my_dpi), dpi=my_dpi) # size of windows

    def update_fig(i):
        print(i)
        # create random data
        x = np.random.rand(100)
        y = np.random.rand(100)
        z = np.random.rand(100)

        fig.clf() # clear figure
        ax=plt.gca()
        ax.tick_params(labelleft=False,labelbottom=False) # remove label
        plt.scatter(x,y,s=1,c=z,vmin=-1,vmax=1,cmap=cm.get_cmap('twilight'),animated=True) # plot something
        plt.tick_params(labelsize=text_fontsize) # axis label size
        plt.xlabel('x',{"fontsize":text_fontsize}) # axis label size
        plt.ylabel('y',{"fontsize":text_fontsize}) # axis label size
        return [fig]
    anim = manimation.FuncAnimation(fig, update_fig,frames=np.arange(begin,end,1), blit=True)
    anim.save("./test.mp4")
    plt.close('all')


if __name__ == '__main__':
    animation()