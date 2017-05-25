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
    print 'Generating J3 function plot with gso_200.out file...'

    x, y = read_coordinates('gso_j3_200.out')

    # First plot    
    fig = plt.figure(figsize=plt.figaspect(0.4))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    X = np.arange(-10, 10, 0.15)
    Y = np.arange(-10, 10, 0.15)
    X, Y = np.meshgrid(X, Y)
    Z = (X*X + Y*Y)**0.25 * (np.sin(50.0*(X*X + Y*Y)**0.1)**2 + 1.0)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=True)
    
    # Legend
    fig.colorbar(surf, shrink=0.5, aspect=10)
    
    # Second plot
    ax = fig.add_subplot(1, 2, 2)
    surf = ax.contourf(X, Y, Z, cmap=cm.coolwarm)
    ax.plot(x, y, '+')

    plt.savefig('gso_j3.png')
    
    print 'gso_j3.png generated.'
