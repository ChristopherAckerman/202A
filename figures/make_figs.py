import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def envelope_theorem(filename='envelope_theorem'):
    x = np.linspace(0, 1, 1000)
    w = -(x - .5)**2 + 1
    v = (-(x - .5)**4) + 1
    fig, ax = plt.subplots()
    ax.plot(x, v, label=r'$v(k)$')
    ax.plot(x, w, label=r'$w(k)$')
    ax.vlines(0.5, ymin=0, ymax=1, linestyle='dashed', color='gray')
    ax.annotate(r'$x_0$', (0.52, 0.03), fontsize=22)
    plt.ylim(0, 1.25)
    plt.legend(loc='lower right', fontsize=16, frameon=False)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.savefig(f'{filename}.pdf')


if __name__ == '__main__':
    envelope_theorem()
