import os
import sys
import logging
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
# -----------------------------------------------------------


class Deproj():
    """
    docstring for Deproj
    """
    def __init__(self, qfile, pixscale = 0.01226, plane = 'disk', nr = 10, nt = 36):
        """
        Class to deproject observations of a debris disk
        """
        self._pixscale = pixscale
        self._nr, self._nt = nr, nt
        self._amin, self._amax, self._incl, self._pa, self._e, self._omega = 0.5, 1.5, 45. * np.pi/180., 0., 0., 0.
        """
        Configure the logging module
        """
        logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s',
            handlers=[
                logging.StreamHandler()
            ])
        self._nx, self._qphi = self._read_fits(qfile)
        """
        Get the coordinates of the center, different if the number
        of pixels is odd or even
        """
        if self._nx%2 == 0:
            self._cx = self._nx //2 - 0.5
        else:
            self._cx = self._nx //2
        self._Xin, self._Yin = (np.mgrid[0:self._nx, 0:self._nx] - self._cx) * self._pixscale
        self._xm, self._ym = np.zeros((self._nx, self._nx)), np.zeros((self._nx, self._nx))
        self._distance, self._azimuth = np.zeros((self._nx, self._nx)), np.zeros((self._nx, self._nx))
        self._xlim = self._cx * self._pixscale

    def go(self, **kwargs):
        """
        Main method to compute the things
        """
        self._check_parameters(kwargs)
        plot = True
        if 'plot' in kwargs:
            plot = kwargs['plot']
        self._xm = ((np.cos(self._pa) * self._Xin + np.sin(self._pa) * self._Yin))
        self._ym = ((np.sin(self._pa) * self._Xin - np.cos(self._pa) * self._Yin)) / np.cos(self._incl)
        self._distance = np.sqrt(self._xm**2. + self._ym**2.)
        self._azimuth = (np.arctan2(self._ym,self._xm) + np.pi/2.) % (2. * np.pi) - np.pi
        # self._azimuth = np.arctan2(self._ym, self._xm)

        rad = np.linspace(self._amin, self._amax, num = self._nr+1)
        theta = np.linspace(-np.pi, np.pi, num = self._nt+1)
        density = np.zeros((self._nr, self._nt))
        for ir in range(self._nr):
            for it in range(self._nt):
                sel = ((self._distance >= rad[ir]) & (self._distance < rad[ir+1]) & (self._azimuth >= theta[it]) & (self._azimuth < theta[it+1]))
                density[ir, it] = np.median(self._qphi[sel])

        if plot:
            """
            Set the cuts
            """
            dmin = np.nanpercentile(density, 1.)
            dmax = np.nanpercentile(density, 99.5)
            qmin = np.nanpercentile(self._qphi, 1.)
            qmax = np.nanpercentile(self._qphi, 99.5)
            """
            Make the plot
            """
            fig = plt.figure(figsize=(15., 5.))
            ax1 = fig.add_subplot(131)
            ax1.imshow(self._qphi, origin='lower', extent = [self._xlim, -self._xlim, -self._xlim, self._xlim], vmin = qmin, vmax = qmax)
            ax1.set_xlabel('$\\alpha$ [$^{\prime\prime}$]')
            ax1.set_ylabel('$\delta$ [$^{\prime\prime}$]')
            ax1.set_xlim(self._amax, -self._amax)
            ax1.set_ylim(-self._amax, self._amax)

            ax2 = fig.add_subplot(132)
            ax2.grid(False)
            ax2.pcolormesh(theta[:-1] * 180. / np.pi, rad[:-1], density, vmin = dmin, vmax = dmax)
            ax2.set_xlabel('$\\theta$ [$^\circ$]')
            ax2.set_ylabel('r [$^{\prime\prime}$]')

            ax3 = fig.add_subplot(133, polar=True)
            ax3.grid(False)
            ax3.set_theta_zero_location('N')
            ax3.pcolormesh(theta[:-1], rad[:-1], density, vmin = dmin, vmax = dmax)
            ax3.set_yticklabels([])
            plt.tight_layout()
            plt.show()

    def debug(self, **kwargs):
        """
        debugging and plot things
        """
        self.go(plot = False, **kwargs)
        r = np.linspace(self._amin, self._amax, num = self._nr)

        fig = plt.figure(figsize=(7,7))
        ax1 = fig.add_axes([0.14, 0.14, 0.8, 0.8])
        ax1.set_aspect('equal')
        ax1.plot(0., 0., marker = '+', color = 'r', ms = 12, mew = 2., zorder = 3)
        cb1 = plt.imshow(self._azimuth, origin = 'lower', extent = [self._xlim, -self._xlim, -self._xlim, self._xlim])
        # plt.contour(self._azimuth, origin = 'lower', extent = [self._xlim, -self._xlim, -self._xlim, self._xlim], levels = self._nt, colors = 'w')
        # cb1 = plt.imshow(self._distance, origin = 'lower', extent = [self._xlim, -self._xlim, -self._xlim, self._xlim])
        for i in range(len(r)):
            xn, yn, theta = self._get_xy(r[i], np.linspace(-np.pi, np.pi, num = self._nt, endpoint = True))
            # ax1.plot(xn, yn, marker = 'o', mec = 'w', mfc = 'w', ms = 4, ls = ' ')
            cb2 = ax1.scatter(xn, yn, c = theta, s = 6)
        # ax1.set_xlabel()
        # ax1.set_ylabel()
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cb1, cax=cax)
        plt.show()

    def _get_xy(self, a, nu):
        """
        All six orbital parameters
        """
        r = a * (1.0 - self._e**2.) / (1.0 + self._e * np.cos(nu))
        x = r * np.cos(nu)
        y = r * np.sin(nu)

        # Now, we need to apply some projection and rotation
        yn = x*(np.cos(self._pa)*np.cos(self._omega) - np.sin(self._omega)*np.cos(self._incl)*np.sin(self._pa)) -\
             y*(np.sin(self._omega)*np.cos(self._pa) +np.cos(self._omega)*np.cos(self._incl)*np.sin(self._pa))
        xn = x*(np.sin(self._pa)*np.cos(self._omega) + np.sin(self._omega)*np.cos(self._incl)*np.cos(self._pa)) +\
             y*(np.cos(self._omega)*np.cos(self._incl)*np.cos(self._pa) - np.sin(self._omega)*np.sin(self._pa))
        # The minus in front of xn is to flip the x axis and put the East on the left
        # The minus in front of nu is to have positive angles on the front side of the disk
        return -xn, yn, -nu

    def _read_fits(self, qfile):
        """
        Method to read the fits file
        """
        if os.path.isfile(qfile):
            hdu = fits.open(qfile)
            data = hdu[0].data
            if len(np.shape(data)) != 2:
                self._error_msg('Problem with the shape of the data: {}'.format(np.shape(data)))
            else:
                if np.shape(data)[0] != np.shape(data)[1]:
                    self._error_msg('Not a square image, not sure if that would work ...')
                else:
                    return np.shape(data)[0], data
        else:
            self._error_msg('Could not find the file')

    def _error_msg(self, message):
        """
        Method to print some error message
        """
        logging.error(message)
        sys.exit(0)

    def _check_parameters(self, kwargs):
        """
        Check parameters that are being passed
        """
        if 'amin' in kwargs:
            self._amin = kwargs['amin']
        if 'amax' in kwargs:
            self._amax = kwargs['amax']
        if 'incl' in kwargs:
            self._incl = kwargs['incl'] * np.pi / 180.
        if 'pa' in kwargs:
            self._pa = -kwargs['pa'] * np.pi / 180.
        if 'e' in kwargs:
            self._e = kwargs['e']
        if 'omega' in kwargs:
            self._omega = kwargs['omega'] * np.pi / 180.

    @property
    def amin(self):
        return self._amin

    @amin.setter
    def amin(self, amin):
        self._amin = amin

    @property
    def amax(self):
        return self._amax

    @amax.setter
    def amax(self, amax):
        self._amax = amax

    @property
    def incl(self):
        return self._incl * 180. / np.pi

    @incl.setter
    def incl(self, incl):
        self._incl = incl * np.pi / 180.

    @property
    def pa(self):
        return self._pa * 180. / np.pi

    @pa.setter
    def pa(self, pa):
        self._pa = pa * np.pi / 180.

    @property
    def e(self):
        return self._e

    @e.setter
    def e(self, e):
        self._e = e

    @property
    def omega(self):
        return self._omega * 180. / np.pi

    @omega.setter
    def omega(self, omega):
        self._omega = omega * np.pi / 180.

#     plt.figure(figsize = (7,7))
#     cb = plt.imshow(distance, origin = 'lower', extent = [xlim, -xlim, -xlim, xlim]) # The axis should now be in arcsec as well
#     plt.contour(distance, origin = 'lower', extent = [xlim, -xlim, -xlim, xlim], colors = 'w', linestyles = '--', levels = [0.5, 1.0, 1.5])
#     plt.plot(xn, yn, color='r')
#     plt.colorbar(cb)
#     plt.show()

"""
Plot the points
"""

if __name__ == '__main__':
    # test = Deproj('test/HD129590_Qphi_300.fits', nr = 11, nt = 30)
    # test.go(amin = 0.15, amax = 0.8, incl = 82.00, pa = -60.61)

    test = Deproj('data_example/HR4796_Qphi_400.fits', nr = 30, nt = 60, pixscale = 0.0072)
    test.go(amin = 0.6, amax = 1.4, incl = 77.72, pa = -151.59)

    # test = Deproj('test/HD121617_Qphi_500.fits', nr = 40, nt = 60)
    # test.go(amin = 0.3, amax = 1.2, incl = 44.6, pa = -118.78)
    # # test.debug(amin = 0.5, amax = 1.5, incl = 80.6, pa = 10., omega = 0.)
