"""
This is the core of the LAURA package and comprehend all subpackages necessary to
do the analysis
Bruno Lustosa de Moura and Leandro de Almeida 2017.
LAURA complete and compact version
"""
from __future__ import print_function, division
import msvcrt
import scipy.fftpack
import numpy as np
import matplotlib.pyplot as plt
import kplr
from time import time as now
import os, sys, errno, math
from numpy.polynomial import Chebyshev as T
from astropy.convolution import Gaussian1DKernel, convolve, Box1DKernel
from astropy.stats import LombScargle
from PyAstronomy.pyTiming import pyPDM
import lmfit
from matplotlib.widgets import Cursor
from astropy.modeling.models import Lorentz1D
import peakutils
from lmfit.models import GaussianModel, LorentzianModel, ExponentialModel
from lmfit import minimize, conf_interval
from matplotlib.widgets import Button
import copy
import scipy as scipy
import pylab as pl
import numpy.random as npr
import random as rd
from time import gmtime, strftime
import glob as glob
from scipy import stats
import logging
plt.style.use('classic')


#--------------------------------------------------------------------------------------------------------------------#
def make_the_timeseries(time, flux):
	assert len(time) == len(flux)
	dt = np.diff(time)
	#print(dt)
	nyquist = (1 / (2 * np.median(dt)))*1e6
	print('Nyquist frequency: %s µHz' % str(nyquist))
	return nyquist
def make_the_power_spectrum(time, flux, oversample):
    clear = lambda:os.system('clear')
    print('Calculating power spectrum')
    freq, power, alpha, beta = power_spectrum(time, flux-1, oversample=oversample)#,memory_use=5000000*10)
    return freq*1e6, power

def cleaning(time,flux,path,limpeza,kic):
    #esta funcao subistitui a noiseremoval
    """
    This function cleans the light curve using the DV
    FLux needs to be normalized at 1 here
    """
    delta_mag = 2*(np.sqrt(2)*np.sqrt(np.mean(np.square(flux-1))))*limpeza
    limSup = 1+delta_mag
    limInf = 1-delta_mag
    xremoved, yremoved = [],[]
    for i in range(len(flux)-1):
        if flux[i]>limSup or flux[i]<limInf:
            yremoved = np.append(yremoved,flux[i])
            xremoved = np.append(xremoved,time[i])
            flux[i]=flux[i-1]
    plt.title('KIC '+np.str(kic))
    plt.xlabel(r'Time [Ms]')
    plt.ylabel(r'Amplitude')
    plt.xlim([time[0], time[-1]])
    plt.subplots_adjust(left=0.12, right=0.95, bottom=0.12, top=0.90)
    plt.plot(time, flux, 'k.')
    plt.plot(xremoved, yremoved, 'rx')
    plt.hlines(limSup,time.min(),time.max(),color="green",linestyles="-")
    plt.hlines(limInf,time.min(),time.max(),color="green",linestyles="-")
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(path+str(kic)+'_TS_cleaned.pdf')
    plt.show()
    return flux
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def get_object(Object,ID,Qtin,Qtfin,method,path):

    client = kplr.API(path)
    if Object == "planet":
        box = client.planet(ID)

    if Object == "star":
        box = client.star(ID)

    lcs = box.get_light_curves(short_cadence=True)
    lcs = lcs[Qtin:Qtfin]

    time, flux_a, flux_c, err_a, err_b = [], [], [], [], []
    s = 0
    for lc in lcs:
        with lc.open() as f:
            hdu_data = f[1].data
            time = np.append(time, hdu_data["time"])
            flux_a = np.append(flux_a, hdu_data["pdcsap_flux"])
            err_a  = np.append(err_a,  hdu_data["pdcsap_flux_err"])
            if s == 0:
                cadencia = (max(time)-min(time))/len(time)
                s = 1
            flux_div, err_div = normalize(hdu_data["pdcsap_flux"], hdu_data["pdcsap_flux_err"], method)
            flux_c = np.append(flux_c, flux_div)
            err_b =  np.append(err_b, err_div)
    return time, flux_a, flux_c, err_a, err_b

def normalize(flux,err,method):
    """ Normalizes the input flux by the method of subtraction or division
         if the division method ("div") is chosen, it will normalize to 1, and if chosen
         the subtraction method ("sub") will normalize by 1e6."""
    onde = 1e8# 671688

    if method == "div":
        mediaF = np.nanmedian(flux)
        mediaE = np.nanmedian(err)
        flux = flux/mediaF
        err = err/mediaE - 1

    if method == "sub":
        mediaF = np.nanmedian(flux)
        flux = flux + (onde - mediaF)

    return flux, err

def detrending(time,flux, err, pl):

    fluxTrue = np.isfinite(flux)
    index = np.where(fluxTrue == True)[0]
    flux_nan = flux[index]
    time = time[index]
    err_nan = err[index]
    time_nan = (time - min(time))

    time_nan, flux_nan = zip(*sorted(zip(time_nan, flux_nan)))
    time_nan = np.array(time_nan)
    flux_nan = np.array(flux_nan)


    p = T.fit(time_nan, flux_nan, pl)
    flux_model = p(time_nan)
    flux_detrended = (flux_nan-flux_model)
    flux_detrended = flux_detrended + 1

    return time_nan, flux_nan, flux_model, flux_detrended, err_nan

def reduction (time_nan, flux_nan, err, cadencia, factor=5, norm=1):


    faltante_time = []
    i = 1
    for i in range (len(time_nan)-1):
        teste = time_nan[i]-time_nan[i-1]

        if teste > 0.0405:
            valor = time_nan[i]-time_nan[i-1]
            falta_time = np.linspace(time_nan[i-1],time_nan[i],(valor/cadencia))
            faltante_time = np.append(faltante_time, falta_time)

    faltante_flux = np.zeros(len(faltante_time))+norm
    faltante_err = np.zeros(len(faltante_time))
    faltante_flux = np.random.normal(faltante_flux,np.std(flux_nan)/factor)
    time_nan = np.append(time_nan, faltante_time)
    flux_nan = np.append(flux_nan, faltante_flux)
    err_nan = np.append(err, faltante_err)

    time_nan, flux_nan = zip(*sorted(zip(time_nan, flux_nan)))
    time_nan = np.array(time_nan)
    flux_nan = np.array(flux_nan)

    return time_nan, flux_nan, err_nan


#-------------------------------------------------------------------------------------------------------------------#
def TS_PS(ID,Qtin,Qtfin,limpeza,oversample,path = 'TEMP/'):

    clear = lambda:os.system('clear')
    path = path+str(ID)+"/"
    mkdir_p(path)
    method = "div"
    detrendy = 1
    preencher = 1
    kic = ID
    if os.path.isfile(path+'/TS_'+str(kic)+'.txt') == False:
        print("Downloading the star ",np.str(ID))
        try:
            time, flux_original, flux, err_original, err = get_object("star", ID, Qtin, Qtfin, method,path)

        except:
            print("maybe your computer is not connected to the internet")
            time.sleep(2)
            sys.exit(1)
        time_nan, flux_nan, flux_model, flux_detrended, err = detrending(time,flux, err, detrendy)
        time_done, flux_done, err = reduction(time_nan,flux_detrended,err,preencher,factor=1)
        time, flux = time_done, flux_done
        flux = cleaning(time,flux,path,limpeza,kic)
        time = time*24*3600 #
        #-------------------------------------------------------------------------------------------------------------------#
        nyquist = make_the_timeseries(time, flux)
        t, flu = time,flux
        file = np.c_[t, flu, err]
        np.savetxt(path+'/TS_'+str(kic)+'.txt', file, header='Time Series frequency, power and error')
    else:
        print('You already have this star in your computer.')

    if os.path.isfile(path+'/PS_'+str(kic)+'.txt') == False:
        print('Analyzing data...')
        data = np.loadtxt(path+'/TS_'+str(kic)+'.txt')
        t = data[:,0]
        flu = data[:,1]
        print('Generating the PowerSpectrum oversample')
        f, p = make_the_power_spectrum(t, flu, oversample=oversample)
        data = np.c_[f,p]
        np.savetxt(path+'/PS_'+str(kic)+'.txt', data, header=' Power_Spectrum_oversample frequency and power')
    else:
        print('You already have the PowerSpectrum of this star.')
        print('If you want to do another analysis, delete the PS_.txt')

    print('The end of the first step')
    print('#-###################################################')
    print('#-###################################################')

#-#########################################################################################################################################################################


def FWHM(X,Y):
	"""returns the error of FWHM"""
	half_max = max(Y) / 2.0
	d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))

	left_idx = np.where(d > 0)[0]
	right_idx = np.where(d < 0)[-1]
	return abs(X[right_idx] - X[left_idx])
def LS(kic,kernerforca,Pbeg,Pend,PeriodMod,path='TEMP/'):

    path = path+str(kic)
    data = np.loadtxt(path+'\TS_'+str(kic)+'.txt')
    time = data[:,0]
    flux = data[:,1]
    time = time/(60*60*24)

    g = Gaussian1DKernel(kernerforca, mode='oversample',factor=5000)
    z = convolve(flux, g, boundary='wrap',normalize_kernel=True)

    plt.plot(time, z, '-b')
    plt.plot(time, flux, '-k',alpha=0.5)
    plt.show()
    print("Entering the analysis by the method: ", PeriodMod)

    if PeriodMod =='LS':

        f = np.linspace(1/Pend,1/Pbeg,50000)
        p = LombScargle(time, z).power(f)
        y =p**2/f

        fig = plt.figure(figsize=(8,5), dpi=130)
        plt.plot(1./f, y,'k', lw=0.5)
        sig1 = np.percentile(y, [1.0])
        sig5 = np.percentile(y, [1])
        xlim = (1./f[0],1./f[-1])
        plt.plot(xlim, [sig1, sig1], ':', c='r',alpha=1, lw=0.5, label='FAP 1\%')
        plt.xlabel("Period[days]")
        plt.ylabel("Power")
        i = peakutils.indexes(y, thres=0.99)
        prot = 1./f[i]
        plt.legend(loc='best')
        plt.savefig(path+'/'+str(kic)+'_LS.pdf')
        plt.show()

    if PeriodMod =='PDM':

        x = time
        y = z

        S = pyPDM.Scanner(minVal=1/Pend, maxVal=1/Pbeg, dVal=0.00001, mode="frequency")
        P = pyPDM.PyPDM(x, y)
        f1, t1 = P.pdmEquiBinCover(5, 1, S)
        f2, t2 = P.pdmEquiBin(10, S)
        fig = plt.figure(figsize=(8,5), dpi=130)
        plt.xlabel("Frequency [days$^{-1}$]")
        plt.ylabel("Theta")
        plt.plot(f1, t1, 'b-',lw=0.5)
        plt.plot(f2, t2,lw=0.5,alpha=0.5)

        j = peakutils.indexes((-t2), thres=0.95)
        aa = f2[j]
        bb = t2[j]

        labels = np.array(np.round(1/f2[j],2))
        texts = []
        for x, y, s in zip(aa, bb, labels):
            texts.append(plt.text(x, y, s, fontsize=7))

        plt.legend(["Multiple sequences", "Equidistant bins"])
        plt.grid(alpha=0.5)
        plt.show()

#-########################################################################################################################

def OBS(kic):
    path = 'TEMP/'+str(kic)+'/'

    color = "#ff7f0e"
    c1 = "#00ff00"
    c2 = '#66cccc'
    c3 = '#cc00ff'
    c4 = '#ee0000'
    labels = {"l0" : '$l=0$', "l1" : '$l=1$', "l2" : '$l=2$', "lz": 'Lorentz peak'}

    ############################################################################################################
    def separacao(freq,psd):
        """This function will calculate the distance between two frequencies
         It is good to automate the calculation of large and small separations"""
        global X
        X = [],[]
        fig = plt.figure(figsize=(17,5), dpi=130)
        plt.plot(freq, psd, 'k-', lw=0.5, alpha=0.5)
        plt.title("Select the high frequency region")
        plt.xlim(min(freq),max(freq))
        plt.ylim(min(psd),max(psd))

        def onclick(event):
            global X
            x = event.xdata
            X = np.append(X,x)
            if len(X) == 1:
                plt.plot((x,x),(min(psd),max(psd)),'r--',alpha=0.5,lw=0.8)
                plt.xlim(min(freq),max(freq))
                plt.ylim(min(psd),max(psd))
                plt.draw()
            if len(X) == 2:
                plt.plot((X[-1],X[-1]),(min(psd),max(psd)),'r--',alpha=0.5,lw=0.8)
                plt.fill_betweenx((min(psd),max(psd)),(X[-2],X[-2]),(X[-1],X[-1]),color='red',alpha=0.2)
                plt.xlim(min(freq),max(freq))
                plt.ylim(min(psd),max(psd))
                plt.draw()
            if len(X)>2:
                plt.clf()
                plt.plot(freq, psd, 'k', lw=0.5, alpha=0.5)
                plt.plot((X[-2],X[-2]),(min(psd),max(psd)),'r--',alpha=0.5,lw=0.8)
                plt.plot((X[-1],X[-1]),(min(psd),max(psd)),'r--',alpha=0.5,lw=0.8)
                plt.fill_betweenx((min(psd),max(psd)),(X[-2],X[-2]),(X[-1],X[-1]),color='red',alpha=0.2)
                plt.xlim(min(freq),max(freq))
                plt.ylim(min(psd),max(psd))
                plt.draw()
            print('Ultimo click: x = ', x)

        fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
        plt.clf()
        return abs(X[-2]-X[-1]), X[-2],X[-1]

    def filtragem(x,y,ajuste):

        plt.subplots_adjust(bottom=0.2)
        box_kernel = Box1DKernel(25)
        yha = convolve(y, box_kernel)
        yhat = ((y - yha)**2)
        yhat = (yhat/yhat.mean())
        yoriginal = copy.copy(yhat)
        original, = plt.plot(x, yoriginal, 'b-',alpha=0.3)
        l, = plt.plot(x, yhat, 'k-',alpha=0.8)
        tresh = 0.05
        j = peakutils.indexes(yhat, thres=tresh)
        inf = min(yhat[j])
        linha, = plt.plot((min(x),max(x)),(inf,inf),ls='--', color = 'blue', lw=0.8)

        k = peakutils.indexes(yhat,thres=tresh)
        xx = x[k]
        yy = yhat[k]
        s = Lorentz1D(xx, yy, fwhm=0.025)
        freqs = np.array(s.amplitude)
        psds = np.array(s.x_0)
        modos, = plt.plot(freqs, psds, 'x', color = 'green',alpha=0.4, label=labels['lz'])
        plt.annotate(np.str(tresh),xy=(np.mean(xx), max(yy)/2))

        class Index(object):
            ajuste = 25
            ajuste2 = 0.05
            def next(self,event):
                global yhat
                self.ajuste += 5
                print (self.ajuste)
                box_kernel = Box1DKernel(self.ajuste)
                yha = convolve(y, box_kernel)
                yhat = ((y - yha)**2)
                yhat = (yhat/yhat.mean())
                original.set_ydata(yoriginal)
                l.set_ydata(yhat)
                plt.draw()

            def prev(self,event):
                global yhat
                self.ajuste -= 5
                if self.ajuste <1:
                    self.ajuste = 5
                print(self.ajuste)
                box_kernel = Box1DKernel(self.ajuste)
                yha = convolve(y, box_kernel)
                yhat = ((y - yha)**2)
                yhat = (yhat/yhat.mean())
                original.set_ydata(yoriginal)
                l.set_ydata(yhat)
                plt.draw()

            def up(self,event):
                global inf, limite, yhat
                self.ajuste2 += 0.01
                j = peakutils.indexes(yhat, thres=self.ajuste2)
                inf = min(yhat[j])
                infinito = (inf,inf)
                k = peakutils.indexes(yhat,thres=self.ajuste2)
                xx = x[k]
                yy = yhat[k]
                s = Lorentz1D(xx, yy, fwhm=0.025)
                freqs = np.array(s.amplitude)
                psds = np.array(s.x_0)
                print('With tresh: ', np.round(self.ajuste2,3), ' Found ', len(psds), ' mods')
                original.set_ydata(yoriginal)
                linha.set_ydata(infinito)
                modos.set_ydata(psds)
                modos.set_xdata(freqs)
                plt.annotate(np.str(self.ajuste2),xy=(np.mean(xx), max(yy)/2))
                limite = self.ajuste2
                plt.draw()

            def down(self,event):
                global inf, limite, yhat
                self.ajuste2 -= 0.01
                if self.ajuste2 <0:
                    self.ajuste2 = 0
                j = peakutils.indexes(yhat, thres=self.ajuste2)
                inf = min(yhat[j])
                infinito = (inf,inf)
                k = peakutils.indexes(yhat,thres=self.ajuste2)
                xx = x[k]
                yy = yhat[k]
                s = Lorentz1D(xx, yy, fwhm=0.025)
                freqs = np.array(s.amplitude)
                psds = np.array(s.x_0)
                print('With tresh ', np.round(self.ajuste2,3), ' found ', len(psds), ' mods')
                original.set_ydata(yoriginal)
                linha.set_ydata(infinito)
                modos.set_ydata(psds)
                modos.set_xdata(freqs)
                limite = self.ajuste2
                plt.draw()


        callback = Index()
        axprev = plt.axes([0.7, 0.05, 0.1, 0.075])
        axnext = plt.axes([0.81, 0.05, 0.1, 0.075])
        axup = plt.axes([0.2, 0.05, 0.1, 0.075])
        axdown = plt.axes([0.31, 0.05, 0.1, 0.075])
        bnext = Button(axnext, 'more')
        bnext.on_clicked(callback.next)
        bprev = Button(axprev, 'less')
        bprev.on_clicked(callback.prev)

        bup = Button(axup, 'up')
        bup.on_clicked(callback.up)

        down = Button(axdown, 'down')
        down.on_clicked(callback.down)
        plt.show()

        return x, yhat, limite



    f = np.loadtxt(str(path)+'PS_'+str(kic)+'.txt')

    freq = f[:,0]
    psd  = f[:,1]

    diff, primeiro, segundo = separacao(freq,psd)

    l = (primeiro<freq)&(freq<segundo)
    x, y = freq[l], psd[l]
    mod = GaussianModel(prefix='gauss_')
    pars = mod.guess(y, x=x)
    out  = mod.fit(y, pars, x=x)
    data = out.fit_report()
    print(data)
    result = out.minimize(method='leastsq')
    ci = conf_interval(out, result, p_names=['gauss_center','gauss_amplitude'], sigmas=(1, 2, 3, 5))
    report = lmfit.printfuncs.report_ci(ci)
    print(report)
    numax = np.array(out.params['gauss_center'])

    plt.annotate('', xy=(primeiro, 0.8*max(psd)), xytext=(primeiro, max(psd)), arrowprops=dict(shrink=0.025, alpha=0.8, fc=c1, ec='k', headwidth=5))
    plt.annotate('', xy=(segundo, 0.8*max(psd)), xytext=(segundo, max(psd)), arrowprops=dict(shrink=0.025, alpha=0.8, fc=c1, ec='k', headwidth=5))
    plt.annotate(r'$\nu_{max}$', xy=(numax, 0.7*max(psd)), xytext=(numax, max(psd)), arrowprops=dict(alpha=0.5, fc='r', ec='r', headwidth=2.4, width=2))
    plt.grid(alpha=0.4)
    plt.ylabel('PSD[Amplitude Units]')
    plt.xlabel('Frequency [$\mu$Hz]')
    plt.title('KIC {}'.format(kic))
    plt.plot(freq, psd, 'k-', lw=0.5, alpha=0.5)
    plt.tight_layout()
    plt.show()
    plt.clf()
    plt.close()

    deltanu = 135*((numax/3050)**0.8)
    np.savetxt(str(path)+'Data_asteroseismic', np.c_[deltanu, numax], header='DeltaNu and Numax')
    print('DeltaNu: {}\nNumax: {}'.format(deltanu, numax))

    j = peakutils.indexes(y, thres=0.015)
    inf = min(y[j])
    plt.axhline(inf,  xmin=0.045, xmax=0.95, ls='--', color='blue', lw=0.8)

    x, yhat, tresh = filtragem(x,y,25)
    print (tresh)
    j = peakutils.indexes(yhat, thres=tresh)
    inf = min(yhat[j])
    plt.axhline(inf,  xmin=0.045, xmax=0.95, ls='--', color='blue', lw=0.8)
    plt.plot(freq[l], yhat, 'k-', lw=0.5)
    k = peakutils.indexes(yhat,thres=tresh)
    xx = x[k]
    yy = yhat[k]
    s = Lorentz1D(xx, yy, fwhm=0.025)
    freqs = np.array(s.amplitude)
    psds = np.array(s.x_0)
    erro = s.fwhm*psds
    dados_lorenz = np.c_[freqs, psds, erro]
    Freq_Region = np.c_[freq[l], yhat]
    np.savetxt(str(path)+str(kic)+'_data_modos_obs.txt', dados_lorenz, header='all peaks (frequencies) of modes detected by lorentz profile, power, error')
    np.savetxt(str(path)+'High_Freq_Region_'+str(kic)+'.txt', Freq_Region, header='High-frequency region with filter Box1DKernel')
    for i in range(0, len(freqs)):
        plt.plot([freqs[i]], [psds[i]], 'x', color = 'green',alpha=0.4, label=labels['lz'])
        labels['lz'] = "_nolegend_"
    plt.tight_layout()
    plt.legend(loc='best',scatterpoints=1)
    plt.savefig(str(path)+'peak_ident_KIC'+str(kic)+'_oversample.png', dpi=170)
    idx = np.argmax(psd[l])
    print('Greater power in the Gaussian region: {}'.format(x[idx]))
    plt.show()

#-################################################################################################################################################

def MOD(kic,largeSep,smallSep,tamanho,Radial_Orders=5):
    """
    Input the star kic
    The biggest separation
    The small separation
    The size of the region comprising each mode (set to 1.5)
    Radial order
    """

    c1 = "#00ff00"
    c2 = '#66cccc'

    path = 'TEMP/'+str(kic)+'/'
    data = np.loadtxt(str(path)+'High_Freq_Region_'+str(kic)+'.txt')
    f, p = data[:,0], data[:,1]

    i = np.argmax(p)
    centralRadial=  f[i]

    Patternfreq = universalPattern(largeSep, smallSep,centralRadial, numberOfRadialOrders = Radial_Orders)
    a = Patternfreq

    i = 0
    radial = a[0]
    dipolo = a[1]
    quadru = a[2]
    ordem = a[4]
    np.savetxt(str(path)+str(kic)+'data_modos_model.txt', np.c_[a[0], a[1], a[2], a[4]], header='frequencies for echelle, radial, di, quadru and the order')

    plt.axvline(40, ls= '-',color=c1, label=r'$l=0$', lw=4,alpha=0.2)
    plt.axvline(40, ls= '-',color=c2, label=r'$l=1$',lw=4,alpha=0.2)

    for i in range(len(radial)):
        plt.axvspan(radial[i]-tamanho/2, radial[i]+tamanho/2, ymin=0.0, ymax=2, alpha=0.2, color=c1)
        plt.axvspan(dipolo[i]-tamanho/2, dipolo[i]+tamanho/2, ymin=0.0, ymax=2, alpha=0.2, color=c2)

    plt.legend(loc='best', numpoints=1)
    plt.plot(f, p, '-r',lw=0.5)
    plt.xlim(f.min(), f.max())
    plt.ylabel('PSD[ppm$^2$/$\mu$Hz]')
    plt.xlabel('Frequency [$\mu$Hz]')
    plt.show()

#-################################################################################################################################
#Aplicacao da biblioteca EcheTed para gerar o diagrama Echelle
"""
Aqui eh um exemplo de como utilizar o pocote EcheTed para gerar um diagrama
Echelle.

Entra:
    Frequencias Modeladas
    Frequencias Obsercadas

Retorna:
    Plot do diagrama
"""

"""
Ted Leandro/Bruno Lustosa @UFRN dez 2017
Rotina para gerar os pontos do diagrama Echelle.
Equacoes usadas de acordo com o seminario
https://meetings.aip.de/592-WE-Heraeus-Seminar/cms/wp-content/uploads/2015/06/Miglio-Bad-Honnef.pdf
"""

def find_l (model,obser,L):
    """
    Encontrar os modos l nas frequencias observadas desorganizadas.
    Entra com:
        model: modos modelados organizados por coluna
        obser: modos observados do jeito que sai do laura

    Retorna:
        lista das frequencias observadas do modo especificado de acordo com
        o modelo

    Eh importante que existam mais modos modelados do que observados.
    """


    l_m = model[:,L]
    l_obser = []
    for i in range(len(l_m)):
        if l_m[i]>int(min(obser)) and l_m[i]<np.round(max(obser),0):
            teste = abs(l_m[i]-obser)
            index = np.where(teste == min(teste))[0]
            valor = obser[index]
            l_obser = np.append(l_obser,float(valor))
        else:
            l_obser = np.append(l_obser,0)
    return l_obser


def organizar(model,obser,qmodos):
    """
    Utiliza o find_l para organizar as frequencias observadas nos modos corretos
    de acordo com as frequencias modelizadas
        model: modos modelados organizados por coluna
        obser: modos observados do jeito que sai do laura

    Retorna:
        array dos modos observados organizados por coluna de acordo com os
        modos que existirem no modelo

    Eh importante que existam mais modos modelados do que observados.
    """
    novo = obser[:,0]
    """
    novo = []   #sera a lista na forma de list das frequencias observadas
    ajustado = 0
    # aqui eu pego todas as colunas do obser e coloco uma embaixo da outra para criar
    #uma lista que sera conferida e organizada
    for i in range(len(obser[0])):
        novo = np.append(novo,obser[:,i]) #coloca as colunas uma embaixo da outra
    """

    if qmodos == 1:
        l_0_obser = np.array(find_l(model,novo,0))

        ajustado = (l_0_obser)

    if qmodos == 2:
        l_0_obser = np.array(find_l(model,novo,0))
        l_1_obser = np.array(find_l(model,novo,1))

        ajustado = (l_0_obser,l_1_obser)

    if qmodos == 3:
        l_0_obser = np.array(find_l(model,novo,0))
        l_1_obser = np.array(find_l(model,novo,1))
        l_2_obser = np.array(find_l(model,novo,2))

        ajustado = (l_0_obser,l_1_obser,l_2_obser)

    if qmodos == 4:
        l_0_obser = np.array(find_l(model,novo,0))
        l_1_obser = np.array(find_l(model,novo,1))
        l_2_obser = np.array(find_l(model,novo,2))
        l_3_obser = np.array(find_l(model,novo,3))

        ajustado = (l_0_obser,l_1_obser,l_2_obser,l_3_obser)

    if qmodos > 4 or qmodos < 1:
        raise ValueError("Rapaz, eh serio que quer organizar "+np.str(qmodos)+" modos? Da nao, cara!")
        exit

    return np.array(ajustado).T

def DELTA(freq):
    """
    Essa funcao recebe as frequencias do modo querido, e retorna o Deltav
    mediano entre essas frequencias. Subtrai-se a seguinte frequencia da anterior
    e cria-se uma lista com essas diferencas, depois, eh feita a media entre essas
    diferencas.

    Entra com:
        frequencias de um modo qualquer

    Retorna:
        DeltaV
    """

    freqDev = []

    for i in range(len(freq)-1):
        freqDev = np.append(freqDev,(freq[i+1]-freq[i]))

    return np.median(freqDev)

def echelle(lm,modo,dev=0):
    """
    As equacoes utilizadas seguem o seminario
    https://meetings.aip.de/592-WE-Heraeus-Seminar/cms/wp-content/uploads/2015/06/Miglio-Bad-Honnef.pdf
    Entra com:
        lm:frequencias em formato de lista
        modo:modo dessas frequencias
        dev:Deltav em questao

    Retorna:
        Eixo x e eixo y do diagrama Echelle na forma vmodDel e frequencia.

    Se a opcao dev for deixada sem nenhuma opcao, o DeltaV sera calculado automaticamente.
    Entao se quiser especificar um Deltav, basta colocar o valor apos o modo.
    """

    l = lm[:,modo]

    nan_index = np.where(l != 0)[0]

    K = nan_index[0]

    ll = l[nan_index]

    vmodDelv = []

    y = []

    if dev == 0:
        delta = DELTA(l)

    else:
        delta = dev
    print("O Deltav utilizado eh ", delta)

    for i in range (len(ll)):
        nu = ll[i]
        k = K+i
        yy = nu
        y = np.append(y,yy)
        vmodv = nu - k*delta
        vmodDelv = np.append(vmodDelv,vmodv)
    return vmodDelv, y


def inclination(N1,Y1,N2,Y2,s,alpha):
    """
    Esta funcao plota os tracos entre os modos e suas frequencias
    Escolher para N1, Y1, N2, Y2 os valores de x e y dos modos mais distantes um do outro
    respectivamente.

    s: o quanto se quer extender as linhas alem do grafico
    alpha: a transparencia das linhas

    Necessita de um plt.show() em algum lugar depois, pois, somente plota.
    """

    for i in range(len(N1)):
        Xa = N1[i]*(1+s)/2 + N2[i]*(1-s)/2
        Ya = Y1[i]*(1+s)/2 + Y2[i]*(1-s)/2
        Xb = N2[i]*(1+s)/2 + N1[i]*(1-s)/2
        Yb = Y2[i]*(1+s)/2 + Y1[i]*(1-s)/2

        plt.plot((Xa,Xb),(Ya,Yb),"k--",alpha=alpha)


def ECHELLE(kic,deltav):

    path = 'TEMP/'+str(kic)+'/'
    data_model = np.loadtxt(str(path)+str(kic)+'data_modos_model.txt')
    data_obser = np.loadtxt(str(path)+str(kic)+'_data_modos_obs.txt')


    data_obser = organizar(data_model,data_obser,3)


    NU0, Y0 = echelle(data_model,0,deltav)
    NU1, Y1 = echelle(data_model,1,deltav)
    NU2, Y2 = echelle(data_model,2,deltav)

    NU0o, Y0o = echelle(data_obser,0,deltav)
    NU1o, Y1o = echelle(data_obser,1,deltav)
    NU2o, Y2o = echelle(data_obser,2,deltav)


    plt.figure(figsize=(8,5))

    inclination(NU2,Y2,NU1,Y1,3,0.2)

    plt.plot(NU0,Y0,"ks",markerfacecolor='white',label="l=0 model",markersize=10)
    plt.plot(NU1,Y1,"ko",markerfacecolor='white',label="l=1 model",markersize=10)
    plt.plot(NU1,Y1,"k-",alpha=0.5)
    plt.plot(NU2,Y2,"k^",markerfacecolor='white',label="l=2 model",markersize=10)

    plt.plot(NU0o,Y0o,"ks",markerfacecolor='black',label="l=0 obser",markersize=10)
    plt.plot(NU1o,Y1o,"ko",markerfacecolor='black',label="l=1 obser",markersize=10)
    plt.plot(NU1o,Y1o,"k-",alpha=0.5)
    plt.plot(NU2o,Y2o,"k^",markerfacecolor='black',label="l=2 obser",markersize=10)

    plt.xlabel("Frequency modulo "+r"$\Delta\nu ="+np.str(np.round(deltav,4))+" (\mu Hz)$")
    plt.ylabel("Frequency "+r"$(\mu Hz)$")

    plt.minorticks_on()
    plt.legend()
    plt.savefig(str(path)+str(kic)+'_echelle_diagram.png')
    plt.show()

    #-######################################################################################################################
""" Aqui fazemos uso do modulo ivs.asteroseismology.redgiantfreqs implementado direto no nosso auxiliar
    usamos o purePressureMode, eh usado dentro do universalPattern para conseguir os modos modelados"""

format  = ''
datefmt = ''

logging.basicConfig(level=logging.INFO,format=format,datefmt=datefmt)
logger = logging.getLogger('IVS.RG')


def purePressureMode(nnn,centralRadial, largeSep, smallSep, epsilon, ell):
      """
      Calculates and returns the frequency of a single pure pressure mode for a spherical (\ell=0,1,2, and 3)[Nota Bene: 3 to be implemented].

      smallSep can be set to an abitrary value if not calculating \ell=2 modes.

          @param nnn: radial order of the mode with respect to the radial order of the central radial mode.
          @type nnn: float or integer
          @param largeSep: large Frequency separation in powerspectrum, in muHz
          @type largeSep: float
          @param smallSep: small Frequency separation in powerspectrum, in muHz
          @type smallSep: float
          @param centralRadial: Frequency of the central radial mode, in muHz
          @type centralRadial: float
          @param epsilon: echelle phase shift of the central Radial mode:  epsilon = centralRadial / largeSep - centralRadial // largeSep
          @type epsilon: float
          @param ell:  spherical degree of a mode
          @type ell: float or integer
          @return frequency in muHz
          @rtype float array
      """

      ell = np.float(ell)
      nnn = np.float(nnn)
      smallSep /= largeSep
      dominantRadialOrderPressure = np.floor(centralRadial / largeSep)
      alphaCorr = 0.015 * largeSep**(-0.32)

      if ell == 0.: constantDEll = 0.0
      if ell == 1.: constantDEll = -(ell/2.+0.025)
      if ell == 2.: constantDEll = smallSep
      if ell == 3.: constantDEll = 0.0

      purePressureModeEll = largeSep * ( (dominantRadialOrderPressure + nnn) + epsilon - constantDEll + alphaCorr/2.*(nnn)**2.)

      return purePressureModeEll





def universalPattern(largeSep, smallSep, centralRadial, numberOfRadialOrders = 3):
    """
    Calculates and returns an array containing frequencies the frequencies of the pure pressure modes
    of the degrees spherical degrees \ell= 0,1,2 and 3 for a defined range of radial orders.

    When slicing, the index number is equal to the spherical degree \ell, universalPattern[0] = all \ell=0, universalPattern[2] = all \ell=2,
    The number of the radial order is given as output[-1]. If desired, the array can also be printed to the screen.

      @param largeSep: large Frequency separation in powerspectrum, in muHz
      @type largeSep: float
      @param smallSep: small Frequency separation in powerspectrum, in muHz
      @type smallSep: float
      @param centralRadial: Frequency of the central radial mode, in muHz
      @type centralRadial: float
      @param numberOfRadialOrders: for how many radial orders above and below the Central Radial Mode should the frequencies be calculated.
      @type numberOfRadialOrders: float or integer

      @return array with frequencies in muHz, [-1] indicates the radial order.
      @rtype float array
    """


    epsilon = centralRadial / largeSep - centralRadial // largeSep
    radialOrderRange = np.arange(-numberOfRadialOrders,numberOfRadialOrders+1,dtype='float')
    dominantRadialOrderPressure = np.floor(centralRadial/largeSep)
    logger.info('dominantRadialOrderPressure {0}'.format(dominantRadialOrderPressure))
    univPatternPerOrder = []

    for nnn in radialOrderRange:
      radial =     purePressureMode(nnn,centralRadial, largeSep, smallSep, epsilon, 0.)
      dipole =     purePressureMode(nnn,centralRadial, largeSep, smallSep, epsilon, 1.)
      quadrupole = purePressureMode(nnn,centralRadial, largeSep, smallSep, epsilon, 2.)
      septupole =  purePressureMode(nnn,centralRadial, largeSep, smallSep, epsilon, 3.)

      univPatternPerOrder.append([radial,dipole,quadrupole,septupole,nnn+dominantRadialOrderPressure])

    univPatternPerOrder = np.array(univPatternPerOrder, dtype='float')


    logger.info('radial\t\tdipole\t\tquadrupole\tseptupole\tRadial Order (Pressure)')
    logger.info('--------------------------------------------------------------------------------------')
    for modeValue in univPatternPerOrder: logger.info('{0:.3f}\t\t{1:.3f}\t\t{2:.3f}\t\t{3:.3f}\t\t{4}'.format(modeValue[0],modeValue[1],modeValue[2],modeValue[3],modeValue[4]))
    logger.info('--------------------------------------------------------------------------------------')

    return univPatternPerOrder.T

def FIT_PS(kic):
    """Essa rotina fita curvas Lorenzianas e Gaussianas no Power Spectrum somente com o kic da estrela"""
    path = 'TEMP/'+str(kic)+'/'
    file = np.loadtxt(str(path)+'PS_'+str(kic)+'.txt')
    f = file[:,0]
    p = file[:,1]

    def fit_gauss(freq,psd):
        global X
        X = [],[]
        fig = plt.figure(figsize=(8,5), dpi=130)
        plt.loglog(freq, psd, 'k', lw=0.5, alpha=0.5)
        plt.title("Select the gaussian fit region")
        plt.ylabel('PSD[ppm$^2$/$\mu$Hz]')
        plt.xlabel('Frequency [$\mu$Hz]')
        plt.xlim(min(freq),max(freq))
        plt.ylim(min(psd),max(psd))
        def onclick(event):
            global X
            x = event.xdata
            X = np.append(X,x)
            if len(X) == 1:
                plt.plot((x,x),(min(psd),max(psd)),'r--',alpha=0.5,lw=0.8)
                plt.xlim(min(freq),max(freq))
                plt.ylim(min(psd),max(psd))
                plt.draw()
            if len(X) == 2:
                plt.fill_betweenx((min(psd),max(psd)),(X[-2],X[-2]),(X[-1],X[-1]),color='red',alpha=0.2)
                plt.plot((X[-1],X[-1]),(min(psd),max(psd)),'r--',alpha=0.5,lw=0.8)
                plt.xlim(min(freq),max(freq))
                plt.ylim(min(psd),max(psd))
                plt.draw()
            if len(X)>2:
                plt.clf()
                plt.title("Selecione a regiao do fit gaussiano")
                plt.ylabel('PSD[ppm$^2$/$\mu$Hz]')
                plt.xlabel('Frequency [$\mu$Hz]')
                plt.xlim(min(freq),max(freq))
                plt.ylim(min(psd),max(psd))

                plt.tight_layout()
                plt.loglog(freq, psd, 'k', lw=0.5, alpha=0.5)
                plt.plot((X[-2],X[-2]),(min(psd),max(psd)),'r--',alpha=0.5,lw=0.8)
                plt.plot((X[-1],X[-1]),(min(psd),max(psd)),'r--',alpha=0.5,lw=0.8)
                plt.fill_betweenx((min(psd),max(psd)),(X[-2],X[-2]),(X[-1],X[-1]),color='red',alpha=0.2)
                plt.draw()
            print('Ultimo click: x = ', x)

        fig.canvas.mpl_connect('button_press_event', onclick)

        plt.tight_layout()
        plt.show()
        plt.clf()

        dados = (X[-2],X[-1])
        return min(dados),max(dados)
    def gauss(x,y):
        print('Selecione a regiao do fit gaussiano quantas vezes for necessario')
        l1,l2 = fit_gauss(x,y)
        index1 = np.where(abs(x-l1) == min(abs(x-l1)))[0]
        index1 = np.int(index1)
        index2 = np.where(abs(x-l2) == min(abs(x-l2)))[0]
        index2 = np.int(index2)
        print (index1, index2)
        x1 = x[index1:index2]
        y1 = y[index1:index2]
        mod = GaussianModel()
        pars = mod.guess(y1, x=x1)
        out = mod.fit(y1,pars, x=x1)
        xgauss = copy.copy(x1)
        ygauss = copy.copy(out.best_fit)
        #-#######
        x0 = x[0:index1]
        y0 = np.zeros(len(y[0:index1]))
        x2 = x[index2:]
        y2 = np.zeros(len(y[index2:]))
        x1 = np.append(x0,x1)
        x1 = np.append(x1,x2)
        y1 = np.append(y0,out.best_fit)
        y1 = np.append(y1,y2)
        print(out.fit_report(min_correl=0.25))


        mod = GaussianModel()
        pars = mod.guess(y1, x=x1)
        out = mod.fit(y1,pars, x=x1)

        #-#######

        return x1, out.best_fit,xgauss,ygauss

    y = p*1.0e12
    x = f
    xx, yy, xgauss, ygauss = gauss(x,y)

    fiti = get_fit(x,y, 'Select the first')
    fiti2 = get_fit(x,y, 'Select the second')

    L1 = LorentzianModel(prefix='L1_')
    pars = L1.guess(y, x=x)
    out = L1.fit(y, pars, x=x)
    result = out.minimize('least_squares')
    plt.clf()
    plt.close('all')
    fig = plt.figure(figsize=(8,5), dpi=130)
    plt.style.use('classic')
    plt.loglog(f, y, "#ff7f0e", lw=0.5)
    plt.plot(x, fiti, 'k--',lw=0.5)
    plt.plot(x, fiti2, 'k--',lw=0.5)
    plt.plot((x[0],x[-1]),(y[-1]/2,y[-1]/2), 'b--',lw=0.5)
    plt.plot(xgauss, ygauss, 'k--',lw=0.5)
    plt.plot(x, fiti+fiti2 - ((fiti+fiti2)/2) + yy, 'r-',lw=0.8)
    plt.title('KIC '+ np.str(kic))
    plt.ylabel('PSD[ppm$^2$/$\mu$Hz]')
    plt.xlabel('Frequency [$\mu$Hz]')
    plt.ylim(1.0e-4,max(y))
    plt.xlim(min(f),max(f))
    plt.tight_layout()
    plt.savefig(str(path)+'PS_KIC'+str(kic)+'_fit.png', dpi=270)
    plt.show()


def get_fit(x,y,ordem):
    global X, X2
    plt.close('all')
    def fit_lorenz(freq,psd):
        global X

        fig = plt.figure(figsize=(8,5), dpi=130)
        plt.loglog(freq, psd, 'k', lw=0.5, alpha=0.5)
        plt.title(str(ordem)+ " Lorenzian fit")
        plt.ylabel('PSD[ppm$^2$/$\mu$Hz]')
        plt.xlabel('Frequency [$\mu$Hz]')
        plt.xlim(min(freq),max(freq))
        plt.ylim(min(psd),max(psd))
        plt.minorticks_on()
        plt.tick_params(direction='in',which='major',top=True,right=True,left=True,length=8,width=1,labelsize=15)
        plt.tick_params(direction='in',which='minor',top=True,right=True,length=5,width=1,labelsize=15)
        plt.tight_layout()
        def onclick(event):
            global X

            plt.clf()
            x = event.xdata
            y = event.ydata
            print('waint...')
            L1 = LorentzianModel(prefix='L1_')
            pars = L1.guess(psd, x=freq)
            pars.update( L1.make_params())
            pars['L1_center'].set(x, min=x-10, max=x+10)
            #pars['L1_sigma'].set(y, min=0)
            pars['L1_amplitude'].set(y, min=0)
            out = L1.fit(psd, pars, x=freq)
            X = out.best_fit
            plt.title(str(ordem)+ " Lorenzian fit")
            plt.ylabel('PSD[ppm$^2$/$\mu$Hz]')
            plt.xlabel('Frequency [$\mu$Hz]')
            plt.minorticks_on()
            plt.tick_params(direction='in',which='major',top=True,right=True,left=True,length=8,width=1,labelsize=15)
            plt.tick_params(direction='in',which='minor',top=True,right=True,length=5,width=1,labelsize=15)
            plt.tight_layout()
            plt.loglog(freq, psd, 'k', lw=0.5, alpha=0.5)
            plt.plot(freq,out.best_fit,'k-', lw=0.5, alpha=0.5)
            plt.xlim(min(freq),max(freq))
            plt.ylim(min(psd),max(psd))
            plt.draw()
            print('Done')


        fig.canvas.mpl_connect('button_press_event', onclick)

        plt.minorticks_on()
        plt.tick_params(direction='in',which='major',top=True,right=True,length=8,width=1,labelsize=15)
        plt.tick_params(direction='in',which='minor',top=True,right=True,length=5,width=1,labelsize=15)
        plt.tight_layout()
        plt.show()
        plt.clf()
        return X

    fiti = fit_lorenz(x,y)
    return fiti

def power_spectrum(time, amplitude, weight=None, minfreq=None, maxfreq=None,
                   oversample=None, memory_use=None, freq=None):
    """
    This function returns the power spectrum of the desired star.
    Arguments:
        - 'time': Time in megaseconds from the timeserie analysis.
        - 'amplitude': Photometry data from the timeserie analysis.
        - 'weight': Weights for each point in the time series.
        - 'minfreq': The lower bound for the frequency interval
        - 'maxfreq': The upper bound for the frequency interval
        - 'oversample': The resolution of the power spectrum.
        - 'memory_use': The amount of memory used for this calculation.
        - 'freq': Override minfreq, maxfreq, ... and use these frequencies instead.
    """

    if minfreq is None:
        minfreq = 1 / (time[-1] - time[0])

    if maxfreq is None:
        maxfreq = 1 / (2 * np.median(np.diff(time)))

    if oversample is None:
        oversample = 4

    if memory_use is None:
        memory_use = 500000
    if weight is None:
        weight = np.ones(amplitude.shape)
    else:
        weight = np.asarray(weight)
        assert weight.shape == amplitude.shape

    if freq is None:

        step = 1 / (oversample * (time[-1] - time[0]))
        freq = np.arange(minfreq, maxfreq, step)
    print("Computing power spectrum for frequencies " +
          "[%g, %g] with step %g" %
          (freq.min(), freq.max(), np.median(np.diff(freq))))


    alpha = np.zeros((len(freq),))
    beta = np.zeros((len(freq),))
    nu = 2 * np.pi * freq
    timerStart = now()
    print_every = 75e6 // len(time)
    chunksize = memory_use // len(time)
    chunksize = max(chunksize, 1)
    print_every = (print_every // chunksize) * chunksize

    for i in range(0, len(nu), chunksize):
        j = min(i + chunksize, len(nu))
        rows = j - i
        if i % print_every == 0:
            elapsedTime = now() - timerStart
            if i == 0:
                totalTime = 0.004 * len(nu)
            else:
                totalTime = (elapsedTime / i) * len(nu)

            print("Progress: %.2f%% (%d of %d)  "
                  "Elapsed: %.2f s  Total: %.2f s"
                  % (np.divide(100.0*i, len(nu)), i, len(nu),
                     elapsedTime, totalTime))

        """
        The outer product is calculated. This way, the product between
        time and ang. freq. will be calculated elementwise; one column
        per frequency. This is done in order to save computing time.
        """
        nutime = np.outer(time, nu[i:j])

        """
        An array with the measured amplitude is made so it has the same size
        as "nutime", since we want to multiply the two.
        """
        amplituderep = amplitude.reshape(-1, 1)
        weightrep = weight.reshape(-1, 1)
        sin_nutime = np.sin(nutime)
        cos_nutime = np.cos(nutime)

        s = np.sum(weightrep * sin_nutime * amplituderep, axis=0)
        c = np.sum(weightrep * cos_nutime * amplituderep, axis=0)
        ss = np.sum(weightrep * sin_nutime ** 2, axis=0)
        cc = np.sum(weightrep * cos_nutime ** 2, axis=0)
        sc = np.sum(weightrep * sin_nutime * cos_nutime, axis=0)

        alpha[i:j] = ((s * cc) - (c * sc)) / ((ss * cc) - (sc ** 2))
        beta[i:j] = ((c * ss) - (s * sc)) / ((ss * cc) - (sc ** 2))

    alpha = alpha.reshape(-1, 1)
    beta = beta.reshape(-1, 1)
    freq = freq.reshape(-1, 1)
    power = alpha ** 2 + beta ** 2
    elapsedTime = now() - timerStart
    print('Computed power spectrum in %.2f s' % (elapsedTime))
    return (freq, power, alpha, beta)


