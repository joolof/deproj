import os
import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
# -----------------------------------------------------------


class Deproj():
    """
    docstring for Deproj
    """
    def __init__(self, qfile, pixscale = 0.01226, nr = 30, nt = 60):
        """
        Class to deproject observations of a debris disk.

        You should provide:
          + name of a fits file (mandatory)
          + pixscale: size of one pixel in arcsec (default 0.01226)
          + nr: number of radial bins (default 10)
          + nt: number of azimuthal bins (default 60)
        """
        self._pixscale = pixscale
        self._nr, self._nt = nr, nt
        self._amin, self._amax, self._incl, self._pa, self._e, self._omega = 0.5, 1.5, 45. * np.pi/180., 0., 0., 0.
        self._nx, self._data = self._read_fits(qfile)
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
        self.density = np.zeros((self._nr, self._nt))
        self._xlim = self._cx * self._pixscale

    def go(self, **kwargs):
        """
        Main method to compute the things

        You should provide the following input parameters:
          + amin: minimum separation (in arcsec)
          + amax: maximum separation (in arcsec)
          + incl: inclination (in degrees)
          + pa: position angle (in degrees)

        Additionally, there are the following optional parameters
          + cmap: a color map (default is viridis)
          + vmin: lower percentile (default 1%)
          + vmax: upper percentile (default 99.5%)
          + xlim: max distance when plotting the observations
          + plot: 'full', 'polar', 'cartesian'
        """
        self._check_parameters(kwargs)
        doplot, plot, cmap = True, 'full', 'viridis'
        vmin, vmax = 1., 99.5
        if 'doplot' in kwargs:
            doplot = kwargs['doplot']
        if 'plot' in kwargs:
            plot = kwargs['plot']
            if (plot != 'full') and (plot != 'polar') and (plot != 'cartesian'):
                self._error_msg('The \'plot\' variable should be either \'full\', \'polar\', or \'cartesian\'.')
        if 'cmap' in kwargs:
            cmap = kwargs['cmap']
        if 'vmin' in kwargs:
            vmin = kwargs['vmin']
        if 'vmax' in kwargs:
            vmax = kwargs['vmax']
        plim = self._xlim
        if 'plim' in kwargs:
            plim = kwargs['plim']
        """
        Do the thing
        """
        self._xm = ((np.cos(self._pa) * self._Xin + np.sin(self._pa) * self._Yin))
        self._ym = ((np.sin(self._pa) * self._Xin - np.cos(self._pa) * self._Yin)) / np.cos(self._incl)
        self._distance = np.sqrt(self._xm**2. + self._ym**2.)
        self._azimuth = (np.arctan2(self._ym,self._xm) + np.pi/2.) % (2. * np.pi) - np.pi
        self._azimuth *= -1 # from clockwise to counter-clockwise because reasons

        rad = np.linspace(self._amin, self._amax, num = self._nr+1)
        theta = np.linspace(-np.pi, np.pi, num = self._nt+1)
        for ir in range(self._nr):
            for it in range(self._nt):
                sel = ((self._distance >= rad[ir]) & (self._distance < rad[ir+1]) & (self._azimuth >= theta[it]) & (self._azimuth < theta[it+1]))
                self.density[ir, it] = np.median(self._data[sel])

        nrad = (rad[:-1] + rad[1:])/2.
        ntheta = (theta[:-1] + theta[1:])/2.

        if doplot:
            """
            Set the cuts
            """
            dmin = np.nanpercentile(self.density, vmin)
            dmax = np.nanpercentile(self.density, vmax)
            qmin = np.nanpercentile(self._data, vmin)
            qmax = np.nanpercentile(self._data, vmax)
            """
            Make the plot
            """
            if plot == 'full':
                fig = plt.figure(figsize=(15., 5.))
                ax1 = fig.add_subplot(131)
                ax1.imshow(self._data, origin='lower', extent = [self._xlim, -self._xlim, -self._xlim, self._xlim], vmin = qmin, vmax = qmax, cmap = cmap)
                ax1.set_xlabel('$\\alpha$ [$^{\prime\prime}$]')
                ax1.set_ylabel('$\delta$ [$^{\prime\prime}$]')
                ax1.set_xlim(plim, -plim)
                ax1.set_ylim(-plim, plim)

                ax2 = fig.add_subplot(132)
                ax2.grid(False)
                ax2.pcolormesh(ntheta * 180./np.pi, nrad, self.density, vmin = dmin, vmax = dmax, cmap = cmap)
                ax2.set_xlabel('$\\theta$ [$^\circ$]')
                ax2.set_ylabel('r [$^{\prime\prime}$]')

                ax3 = fig.add_subplot(133, polar=True)
                ax3.grid(False)
                ax3.set_theta_zero_location('N')
                ax3.set_theta_direction(-1)
                ax3.pcolormesh(ntheta, nrad, self.density, vmin = dmin, vmax = dmax, cmap = cmap)
                ax3.set_xticks(np.linspace(np.pi, -np.pi, 8, endpoint=False))
                ax3.set_thetalim(-np.pi, np.pi)
                ax3.set_yticklabels([])
                plt.tight_layout()
                plt.show()
            elif plot == 'polar':
                fig = plt.figure(figsize=(7,7))
                ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar = True)
                ax1.grid(False)
                ax1.set_theta_zero_location('N')
                ax1.set_theta_direction(-1)
                ax1.pcolormesh(ntheta, nrad, self.density, vmin = dmin, vmax = dmax, cmap = cmap)
                ax1.set_xticks(np.linspace(np.pi, -np.pi, 8, endpoint=False))
                ax1.set_thetalim(-np.pi, np.pi)
                ax1.set_yticklabels([])
                plt.show()
            else:
                fig = plt.figure(figsize=(7,6))
                ax1 = fig.add_axes([0.14, 0.14, 0.8, 0.8])
                ax1.grid(False)
                ax1.pcolormesh(ntheta * 180. / np.pi, nrad, self.density, vmin = dmin, vmax = dmax, cmap = cmap)
                ax1.set_xlabel('$\\theta$ [$^\circ$]')
                ax1.set_ylabel('r [$^{\prime\prime}$]')
                plt.show()

        return ntheta, nrad, self.density

    def debug(self, **kwargs):
        """
        Method to plot the midplane distance and azimuthal angle.

        You should provide the following input parameters:
          + amin: minimum separation (in arcsec)
          + amax: maximum separation (in arcsec)
          + incl: inclination (in degrees)
          + pa: position angle (in degrees)
        """
        self.go(doplot = False, **kwargs) # Just to get some of the variables

        fig = plt.figure(figsize=(10., 4.5))
        ax1 = fig.add_subplot(121)
        ax1.set_aspect('equal')
        ax1.plot(0., 0., marker = '+', color = 'w', ms = 8, mew = 2., zorder = 3)
        cb1 = ax1.imshow(self._distance, origin = 'lower', extent = [self._xlim, -self._xlim, -self._xlim, self._xlim])
        CS = ax1.contour(self._distance, origin = 'lower', extent = [self._xlim, -self._xlim, -self._xlim, self._xlim], levels=[self._amin, self._amax], colors = 'w', linewidths=1)
        ax1.clabel(CS, inline=True, fontsize=8)
        ax1.set_xlabel('$\\alpha$ [$^{\prime\prime}$]')
        ax1.set_ylabel('$\delta$ [$^{\prime\prime}$]')
        ax1.set_title('Midplane distance [$^{\prime\prime}$]', fontdict={'fontsize':10})
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cb1, cax=cax)

        ax2 = fig.add_subplot(122)
        ax2.set_aspect('equal')
        ax2.plot(0., 0., marker = '+', color = 'w', ms = 8, mew = 2., zorder = 3)
        cb1 = ax2.imshow(self._azimuth * 180./np.pi, origin = 'lower', extent = [self._xlim, -self._xlim, -self._xlim, self._xlim])
        ax2.set_xlabel('$\\alpha$ [$^{\prime\prime}$]')
        ax2.set_title('Midplane azimuth [$^\circ$]', fontdict={'fontsize':10})
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cb1, cax=cax)
        plt.tight_layout()
        plt.show()

    def _read_fits(self, qfile):
        """
        Method to read the fits file
        """
        if type(qfile) is str:
            if os.path.isfile(qfile):
                hdu = fits.open(qfile)
                data = hdu[0].data
                hdu.close()
                if len(np.shape(data)) != 2:
                    self._error_msg('Problem with the shape of the data: {}'.format(np.shape(data)))
                else:
                    if np.shape(data)[0] != np.shape(data)[1]:
                        self._error_msg('Not a square image, not sure if that would work ...')
                    else:
                        return np.shape(data)[0], data
            else:
                self._error_msg('Could not find the file')
        elif type(qfile) is np.ndarray:
            data = qfile
            if len(np.shape(data)) != 2:
                self._error_msg('Problem with the shape of the data: {}'.format(np.shape(data)))
            else:
                if np.shape(data)[0] != np.shape(data)[1]:
                    self._error_msg('Not a square image, not sure if that would work ...')
                else:
                    return np.shape(data)[0], data
        else:
            self._error_msg('The data should either be the name of a fits file or a 2D array.')

    def _error_msg(self, message):
        """
        Method to print some error message
        """
        print(message)
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
        return -self._pa * 180. / np.pi

    @pa.setter
    def pa(self, pa):
        self._pa = -pa * np.pi / 180.

"""
Plot the points
"""

if __name__ == '__main__':
    # test = Deproj('test/HD129590_Qphi_300.fits', nr = 11, nt = 30)
    # test.go(amin = 0.15, amax = 0.8, incl = 82.00, pa = -60.61, plot = 'cartesian')

    # test = Deproj('data_example/HR4796_Qphi_400.fits', nr = 30, nt = 60, pixscale = 0.0072)
    # test.go(amin = 0.7, amax = 1.3, incl = 77.72, pa = -151.59, plot = 'polar')
    # test.debug(amin = 0.7, amax = 1.3, incl = 77.72, pa = -151.59)

    test = Deproj('test/HD121617_Qphi_500.fits', nr = 30, nt = 60)
    test.go(amin = 0.4, amax = 1.2, incl = 44.6, pa = -118.78, plot = 'full', plim = 1.2)
    # test.debug(amin = 0.4, amax = 1.2, incl = 44.6, pa = -118.78)
