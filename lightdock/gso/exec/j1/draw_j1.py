import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm


def read_coordinates(file_name):
    """Reads a GSO output population file and returns x and y coordinates."""
    f = open(file_name)
    x = []
    y = []
    for line in f:
        if not line.startswith('#'):
            i = line.index('(')
            j = line.index(')')
            raw = line[i+1:j]
            values = raw.split(',')
            x.append(float(values[0]))
            y.append(float(values[1]))
    return x, y


if __name__ == "__main__":
    print 'Generating J1 function plot with gso_200.out file...'

    x, y = read_coordinates('gso_j1_200.out')

    # First plot    
    fig = plt.figure(figsize=plt.figaspect(0.4))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    X = np.arange(-3, 3, 0.15)
    Y = np.arange(-3, 3, 0.15)
    X, Y = np.meshgrid(X, Y)
    Z = 3.0 * (1-X)*(1-X) * np.exp(-(X*X + (Y+1)*(Y+1))) \
        - 10.0 * (X/5.0 - (X**3) - (Y**5)) * np.exp(-(X*X + Y*Y)) - 1/3.0 * np.exp(-((X+1)*(X+1) + Y*Y))
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
    
    # Legend
    fig.colorbar(surf, shrink=0.5, aspect=10)
    
    # Second plot
    ax = fig.add_subplot(1, 2, 2)
    surf = ax.contourf(X, Y, Z, cmap=cm.coolwarm)
    ax.plot(x, y, '+')

    plt.savefig('gso_j1.png')
    
    print 'gso_j1.png generated.'
