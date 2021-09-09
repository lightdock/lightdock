"""Module to calculate the best fitting ellipsoid given a set of 3D points.

Adapted from https://github.com/minillinim/ellipsoid/blob/master/ellipsoid.py
"""

import numpy as np
from numpy import linalg
from lightdock.error.lightdock_errors import MinimumVolumeEllipsoidError


class MinimumVolumeEllipsoid(object):
    def __init__(self, points, precision=0.01):
        self.center = None
        self.radii = None
        self.rotation = None
        self.poles = None
        self.axes = None
        self._points = points
        self._precision = precision
        try:
            self._get_min_vol_ellipsoid()
        except Exception as e:
            raise MinimumVolumeEllipsoidError(
                "Can not build minimum volume ellipsoid. Reason: %s" % str(e)
            )

    def _get_min_vol_ellipsoid(self):
        """Finds the minimum volume ellipsoid which contains the set of given points"""
        (N, d) = np.shape(self._points)
        d = float(d)

        # Q will be our working array
        Q = np.vstack([np.copy(self._points.T), np.ones(N)])
        QT = Q.T

        # initializations
        err = 1.0 + self._precision
        u = (1.0 / N) * np.ones(N)

        # Khachiyan Algorithm
        while err > self._precision:
            V = np.dot(Q, np.dot(np.diag(u), QT))
            M = np.diag(
                np.dot(QT, np.dot(linalg.inv(V), Q))
            )  # M the diagonal vector of an NxN matrix
            j = np.argmax(M)
            maximum = M[j]
            step_size = (maximum - d - 1.0) / ((d + 1.0) * (maximum - 1.0))
            new_u = (1.0 - step_size) * u
            new_u[j] += step_size
            err = np.linalg.norm(new_u - u)
            u = new_u

        # center of the ellipsoid
        self.center = np.dot(self._points.T, u)

        # the A matrix for the ellipsoid
        A = (
            linalg.inv(
                np.dot(self._points.T, np.dot(np.diag(u), self._points))
                - np.array([[a * b for b in self.center] for a in self.center])
            )
            / d
        )

        _, s, self.rotation = linalg.svd(A)
        # s can be zero in degenerated cases causing RuntimeWarning divide by zero error
        np.seterr(all="raise")
        self.radii = 1.0 / np.sqrt(s)

        # Calculate poles
        axes = np.array(
            [
                [self.radii[0], 0.0, 0.0],
                [0.0, self.radii[1], 0.0],
                [0.0, 0.0, self.radii[2]],
            ]
        )

        for i in range(len(axes)):
            axes[i] = np.dot(axes[i], self.rotation)
        self.axes = axes

        self.poles = []
        for p in axes:
            x1 = -p[0] + self.center[0]
            x2 = p[0] + self.center[0]
            y1 = -p[1] + self.center[1]
            y2 = p[1] + self.center[1]
            z1 = -p[2] + self.center[2]
            z2 = p[2] + self.center[2]
            self.poles.append([x1, y1, z1])
            self.poles.append([x2, y2, z2])
