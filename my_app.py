from __future__ import annotations
from PyQt5 import QtWidgets, QtGui, uic
from PyQt5 import QtCore
#from pyqtgraph import PlotWidget, plot
#import pyqtgraph as pg
import sys  # We need sys so that we can pass argv to QApplication
import os
import numpy as np
from typing import List, Dict, Tuple, Optional, Type, Any, Union

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # Load the UI Page
        uic.loadUi('mainwindow.ui', self)

        palr = self.label_camera1L.palette()
        palg = self.label_camera1R.palette()
        palr.setColor(QtGui.QPalette.WindowText, QtGui.QColor("red"))
        palg.setColor(QtGui.QPalette.WindowText, QtGui.QColor("green"))
        self.label_camera1L.setPalette(palr)
        self.label_camera1R.setPalette(palr)
        self.label_camera2L.setPalette(palg)
        self.label_camera2R.setPalette(palg)

        # define default values for objectives and update ui elements
        objmag = 20
        objna = 0.95
        addmag = 1.0

        # define other default values and update ui elements
        emwl = 520
        phf = 50
        sampling = 2.0

        # store them
        self.objmag1.setValue(objmag)
        self.objmag2.setValue(objmag)
        self.objna1.setValue(objna)
        self.objna2.setValue(objna)
        self.addmag1.setValue(addmag)
        self.addmag2.setValue(addmag)

        # define default values for optics etc.
        self.mic1 = Microscope(name="Mic1",
                               objmag=objmag,
                               objna=objna,
                               addmag=addmag)

        self.mic2 = Microscope(name="Mic2",
                               objmag=objmag,
                               objna=objna,
                               addmag=addmag)

        self.emwl_value = emwl
        self.emwl.setValue(emwl)
        self.phf1_value = phf
        self.phf1.setValue(phf)
        self.sampling_value = sampling
        self.nyq.setValue(sampling)

        # initialize tow cameras with default values
        self.cam1 = Camera(name="CAM1",
                           qe=0.75,
                           pixsize=4.54,
                           binning=1,
                           cameratype="CCD",
                           emgain=1,
                           readout=6.0,
                           dark=0.005,
                           cic=0.0)

        self.cam2 = Camera(name="CAM2",
                           qe=0.65,
                           pixsize=6.45,
                           binning=1,
                           cameratype="CCD",
                           emgain=1,
                           readout=7.0,
                           dark=0.005,
                           cic=0.0)

        # adapt the noise factor and readout noise
        self.cam1 = adapt_readout(self.cam1)
        self.cam1 = adapt_noisefactor(self.cam1)
        self.cam2 = adapt_readout(self.cam2)
        self.cam2 = adapt_noisefactor(self.cam2)

        # update the noise and size factors
        self.noisef1.setText(str(self.cam1.nf))
        self.noisef2.setText(str(self.cam2.nf))
        self.sizef1.setText(str(self.cam1.nf))
        self.sizef2.setText(str(self.cam2.nf))

        # calculate the values for both cameras
        self.cp1, self.cp2 = calc_values(self.cam1, self.cam2, self.mic1, self.mic2,
                                         emwl=self.emwl_value,
                                         phf=self.phf1_value,
                                         sampling=self.sampling_value)

        print("Cameras initialized and values calculated.")

        # update ui
        self.phf2.setText(str(self.cp2["flux"]))

        # update values for the pixel sizes
        self.piximage1.setText(str(self.cp1["piximage"]))
        self.piximage2.setText(str(self.cp2["piximage"]))
        self.pixrequired1.setText(str(self.cp1["req_pixsize"]))
        self.pixrequired2.setText(str(self.cp2["req_pixsize"]))

        # configure plot
        self.MplWidget.canvas.axes.set_title("Camera SNR Plot")
        self.MplWidget.canvas.axes.set_xlabel("Photons / Pixel / Frame")
        self.MplWidget.canvas.axes.set_ylabel("SNR Ratio")
        self.MplWidget.canvas.axes.grid(True)
        self.MplWidget.canvas.axes.set_xlim(0, 200)
        self.MplWidget.canvas.axes.set_ylim(0, 10)

        # plot SNR curves (returns a tuple of line objects, thus the comma)
        self.snr1_curve, = self.MplWidget.canvas.axes.plot(self.cp1["phf"], self.cp1["snr"], "r-", lw=3)
        self.snr2_curve, = self.MplWidget.canvas.axes.plot(self.cp1["phf"], self.cp2["snr"], "g-", lw=3)

        # plot indicator lines (returns a tuple of line objects, thus the comma)
        self.indicator1_line, = self.MplWidget.canvas.axes.plot(self.cp1["phindx"], self.cp1["phindy"], "r-.", lw=2)
        self.indicator2_line, = self.MplWidget.canvas.axes.plot(self.cp2["phindx"], self.cp2["phindy"], "g-.", lw=2)

        # connect scaling values for the plot
        self.xscale_min.valueChanged.connect(self.change_scale)
        self.xscale_max.valueChanged.connect(self.change_scale)
        self.yscale_min.valueChanged.connect(self.change_scale)
        self.yscale_max.valueChanged.connect(self.change_scale)

        # connect binning selectors
        self.bin1.currentIndexChanged.connect(self.change_binning)
        self.bin2.currentIndexChanged.connect(self.change_binning)

        # connect qe values
        self.qe1.valueChanged.connect(self.change_qe)
        self.qe2.valueChanged.connect(self.change_qe)

        # connect pixel size values
        self.pixsize1.valueChanged.connect(self.change_pix)
        self.pixsize2.valueChanged.connect(self.change_pix)

        # connect readout noise values
        self.readnoise1.valueChanged.connect(self.change_readoutnoise)
        self.readnoise2.valueChanged.connect(self.change_readoutnoise)

        # connect camera type values
        self.type1.currentIndexChanged.connect(self.change_type)
        self.type2.currentIndexChanged.connect(self.change_type)

        # connect emgain values
        self.emgain1.valueChanged.connect(self.change_gain)
        self.emgain2.valueChanged.connect(self.change_gain)

        # connect dark current values
        self.dark1.valueChanged.connect(self.change_dark)
        self.dark2.valueChanged.connect(self.change_dark)

        # connect the CIC noise values
        self.cic1.valueChanged.connect(self.change_cic)
        self.cic2.valueChanged.connect(self.change_cic)

        # connect objective magnification values
        self.objmag1.valueChanged.connect(self.change_objmag)
        self.objmag2.valueChanged.connect(self.change_objmag)

        # connect additional magnification values
        self.objna1.valueChanged.connect(self.change_objna)
        self.objna2.valueChanged.connect(self.change_objna)

        # connect additional magnification values
        self.addmag1.valueChanged.connect(self.change_addmag)
        self.addmag2.valueChanged.connect(self.change_addmag)

        # connect sampling value
        self.nyq.valueChanged.connect(self.change_sampling)

        # connect EM-WL  value
        self.emwl.valueChanged.connect(self.change_emwl)


    # modify plot

    def change_scale(self: QtWidgets.QMainWindow) -> None:

        self.MplWidget.canvas.axes.set_xlim(self.xscale_min.value(), self.xscale_max.value())
        self.MplWidget.canvas.axes.set_ylim(self.yscale_min.value(), self.yscale_max.value())

        # update the plot
        self.MplWidget.canvas.draw()

    # change camera parameters

    def change_binning(self: QtWidgets.QMainWindow) -> None:

        # change binning values
        self.cam1.binning = self.bin1.currentIndex() + 1
        self.cam2.binning = self.bin2.currentIndex() + 1

        # update the plot and redraw
        self.update_plot()

    def change_qe(self: QtWidgets.QMainWindow) -> None:

        # change the qe values values
        self.cam1.qe = self.qe1.value()
        self.cam2.qe = self.qe2.value()

        # update the plot and redraw
        self.update_plot()

    def change_pix(self: QtWidgets.QMainWindow) -> None:

        # change the pixel size values
        self.cam1.pixsize = self.pixsize1.value()
        self.cam2.pixsize = self.pixsize2.value()

        # update the plot and redraw
        self.update_plot()

    def change_readoutnoise(self: QtWidgets.QMainWindow) -> None:

        # change the readout noise
        self.cam1.readout = self.readnoise1.value()
        self.cam2.readout = self.readnoise2.value()

        # update the plot and redraw
        self.update_plot()

    def change_dark(self: QtWidgets.QMainWindow) -> None:

        # change the readout noise
        self.cam1.dark = self.dark1.value()
        self.cam2.dark = self.dark2.value()

        # update the plot and redraw
        self.update_plot()

    def change_cic(self: QtWidgets.QMainWindow) -> None:

        # change the readout noise
        self.cam1.cic = self.cic1.value()
        self.cam2.cic = self.cic2.value()

        # update the plot and redraw
        self.update_plot()

    def change_gain(self: QtWidgets.QMainWindow) -> None:
        # change the camera gain
        self.cam1.emgain = self.emgain1.value()
        self.cam2.emgain = self.emgain2.value()

        # update the plot and redraw
        self.update_plot()

    def change_type(self: QtWidgets.QMainWindow) -> None:

        # change the camera type
        self.cam1.cameratype = self.type1.currentText()
        self.cam2.cameratype = self.type2.currentText()

        self.cam1.adapt_noisefactor()
        self.cam1.adapt_readout()
        self.cam2.adapt_noisefactor()
        self.cam2.adapt_readout()

        # update the plot and redraw
        self.update_plot()

    # change optics

    def change_addmag(self: QtWidgets.QMainWindow) -> None:

        # change the readout noise
        self.mic1.addmag = self.addmag1.value()
        self.mic2.addmag = self.addmag2.value()

        # update the plot and redraw
        self.update_plot()

    def change_objmag(self: QtWidgets.QMainWindow) -> None:

        # change the readout noise
        self.mic1.objmag = self.objmag1.value()
        self.mic2.objmag = self.objmag2.value()

        # update the plot and redraw
        self.update_plot()

    def change_objna(self: QtWidgets.QMainWindow) -> None:

        # change the readout noise
        self.mic1.objna = self.objna1.value()
        self.mic2.objna = self.objna2.value()

        # update the plot and redraw
        self.update_plot()

    def change_sampling(self: QtWidgets.QMainWindow) -> None:

        # change the sampling value
        self.sampling_value = self.nyq.value()

        # update the plot and redraw
        self.update_plot()

    def change_emwl(self: QtWidgets.QMainWindow) -> None:

        # change the readout noise
        self.emwl_value = self.emwl.value()

        # update the plot and redraw
        self.update_plot()


    # update the plot

    def update_plot(self):

        # recalculate the values
        self.cp1, self.cp2 = calc_values(self.cam1, self.cam2, self.mic1, self.mic2,
                                         emwl=self.emwl_value,
                                         phf=self.phf1_value,
                                         sampling=self.sampling_value)

        # update the line data for the SNR curves
        self.snr1_curve.set_xdata(self.cp1["phf"])
        self.snr1_curve.set_ydata(self.cp1["snr"])
        self.snr2_curve.set_xdata(self.cp1["phf"])
        self.snr2_curve.set_ydata(self.cp2["snr"])

        # update the line data for the indicator lines
        self.indicator1_line.set_xdata(self.cp1["phindx"])
        self.indicator1_line.set_ydata(self.cp1["phindy"])
        self.indicator2_line.set_xdata(self.cp2["phindx"])
        self.indicator2_line.set_ydata(self.cp2["phindy"])

        # update values for the pixel sizes
        #self.piximage1.setText(str(self.cp1["piximage"]))

        # update the whole plot
        self.MplWidget.canvas.draw()
        self.MplWidget.canvas.flush_events()


class Camera:
    def __init__(self, name: str = "Camera1",
                 qe: float = 0.75,
                 pixsize: float = 4.54,
                 binning: int = 1,
                 cameratype: str = "CCD",
                 emgain: int = 1,
                 readout: float = 1.0,
                 noisefactor: float = 1.0,
                 dark: float = 0.005,
                 cic: float = 0.0) -> None:

        self.qe = qe
        self.pixsize = pixsize
        self.binning = binning
        self.cameratype = cameratype
        self.emgain = emgain
        self.readout = readout
        self.readout_mod = readout
        self.nf = noisefactor
        self.dark = dark
        self.cic = cic


class Microscope:
    def __init__(self, name: str = "Mic1",
                 objmag: float = 20.0,
                 objna: float = 0.95,
                 addmag: float = 1.0) -> None:

        self.mic1 = name
        self.objmag = objmag
        self.objna = objna
        self.addmag = addmag


def calc_values(cam1: type[Camera], cam2: type[Camera], mic1: type[Microscope], mic2: type[Microscope],
                emwl: int = 520,
                phf: int = 50,
                sampling: float = 2.0) -> Tuple[Dict, Dict]:

    cp1 = {}
    cp2 = {}
    cp1["flux"] = phf

    # pixel size in image plane incl. binning
    cp1["piximage"] = cam1.pixsize * cam1.binning / (mic1.objmag * mic1.addmag)
    cp2["piximage"] = cam2.pixsize * cam2.binning / (mic2.objmag * mic2.addmag)

    # required pixel since in image to fulfil Nyquist
    cp1["req_pixsize"] = float(np.round(0.61 * (emwl / 1000) / (sampling * mic1.objna), 3))
    cp2["req_pixsize"] = float(np.round(0.61 * (emwl / 1000) / (sampling * mic2.objna), 3))

    # correction factor for pixel area
    cp2["corrf_pixarea"] = 1.00
    cp2["corrf_pixarea"] = float(np.round((cp2["piximage"] **2) / (cp1["piximage"] **2), 2))

    # adapt the readout noise if camera is an CMOS
    #if cam1.cameratype == "CMOS":
    #    cp1["readout_mod"] = np.sqrt(cam1.binning)
    #else:
    #    cp1["readout_mod"] = cam1.readout
    #if cam2.cameratype == "CMOS":
    #    cp2["readout_mod"] = np.sqrt(cam2.binning)
    #else:
    #    cp2["readout_mod"] = cam2.readout

    # create ph vector containing the number of detected photons and use for both cameras
    cp1["phf"] = np.arange(0, 400, 1, dtype=np.int16)
    #cp2["phf"] = np.arange(0, 400, 1, dtype=np.int16)

    # calculation of SNR including CIC - Clock Induced Charge
    #cp1["snr"] = (cam1.qe * cp1["phf"] / np.sqrt(cam1.nf**2 * (cam1.qe * cp1["phf"] + cam1.dark**2 + cam1.cic**2) + (cp1["readout_mod"]**2 / cam1.emgain**2))).astype(float)
    #cp2["snr"] = (cam2.qe * cp1["phf"] / np.sqrt(cam2.nf**2 * (cam2.qe * cp1["phf"] + cam2.dark**2 + cam2.cic**2) + (cp2["readout_mod"]**2 / cam2.emgain**2))).astype(float)

    cp1["snr"] = (cam1.qe * cp1["phf"] / np.sqrt(cam1.nf**2 * (cam1.qe * cp1["phf"] + cam1.dark**2 + cam1.cic**2) + (cam1.readout_mod**2 / cam1.emgain**2))).astype(float)
    cp2["snr"] = (cam2.qe * cp1["phf"] / np.sqrt(cam2.nf**2 * (cam2.qe * cp1["phf"] + cam2.dark**2 + cam2.cic**2) + (cam2.readout_mod**2 / cam2.emgain**2))).astype(float)

    # calculate values for photon indicators
    cp2["flux"] = (np.round(cp1["flux"] * cp2["corrf_pixarea"], 0)).astype(int)

    #  calculate explicit SNR values
    #cp1["snr_value"] = ((cam1.qe * cp1["flux"]) / np.sqrt(cam1.nf**2 * (cam1.qe * cp1["flux"] + cam1.dark**2 + cam1.cic**2) + (cp1["readout_mod"]**2 / cam1.emgain**2))).astype(float)
    #cp2["snr_value"] = ((cam2.qe * cp2["flux"]) / np.sqrt(cam2.nf**2 * (cam2.qe * cp2["flux"] + cam2.dark**2 + cam2.cic**2) + (cp2["readout_mod"]**2 / cam2.emgain**2))).astype(float)
    cp1["snr_value"] = ((cam1.qe * cp1["flux"]) / np.sqrt(cam1.nf**2 * (cam1.qe * cp1["flux"] + cam1.dark**2 + cam1.cic**2) + (cam1.readout_mod**2 / cam1.emgain**2))).astype(float)
    cp2["snr_value"] = ((cam2.qe * cp2["flux"]) / np.sqrt(cam2.nf**2 * (cam2.qe * cp2["flux"] + cam2.dark**2 + cam2.cic**2) + (cam2.readout_mod**2 / cam2.emgain**2))).astype(float)

    cp1["phindx"] = np.array([cp1["flux"], cp1["flux"], 0])
    cp1["phindy"] = np.array([0, cp1["snr_value"], cp1["snr_value"]])

    cp2["phindx"] = np.array([cp2["flux"], cp2["flux"], 0])
    cp2["phindy"] = np.array([0, cp2["snr_value"], cp2["snr_value"]])

    return cp1, cp2


def adapt_noisefactor(cam: Camera) -> Camera:
    # adjust noise factor due to CCD type
    if cam.cameratype == "CCD":
        # reset noise factor and gain in case of an normal CCD
        cam.nf = 1.0
        cam.emgain = 1
        cam.cic = 0.0

    elif cameratype == "EM-CCD":
        cam.nf = 1.41

    return cam


def adapt_readout(cam: Camera) -> Camera:
    # adapt the readout noise if camera is an CMOS
    if cam.cameratype == "CMOS":
        cam.readout_mod = np.sqrt(cam.binning)

    return cam


def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()