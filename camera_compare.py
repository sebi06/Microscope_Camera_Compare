# -*- coding: utf-8 -*-

#################################################################
# File        : camera_compare.py
# Author      : sebi06
#
#
# This code probably does not reflect the latest new technologies
# of microscope cameras anymore but is hopefully still useful to compare cameras
# and understand the lines of reasoning when choosing the "right" camera
#
# Disclaimer: The code is purely experimental. Feel free to
# use it at your own risk.
#
#################################################################

from __future__ import annotations
from PyQt5 import QtWidgets, QtGui, uic
from PyQt5 import QtCore
import sys
import os
import numpy as np
from typing import List, Dict, Tuple, Optional, Type, Any, Union


class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        # Load the UI Page
        uic.loadUi("mainwindow.ui", self)

        # on eway to modify the color
        palr = self.label_camera1L.palette()
        palg = self.label_camera1R.palette()
        palr.setColor(QtGui.QPalette.WindowText, QtGui.QColor("red"))
        palg.setColor(QtGui.QPalette.WindowText, QtGui.QColor("green"))
        self.label_camera1L.setPalette(palr)
        self.label_camera1R.setPalette(palr)
        self.label_camera2L.setPalette(palg)
        self.label_camera2R.setPalette(palg)

        # another way to modify the color
        self.name1.setStyleSheet("""QLineEdit {color: red }""")
        self.name2.setStyleSheet("""QLineEdit {color: green }""")
        self.phf1.setStyleSheet("""QSpinBox {color: red }""")
        self.phf2.setStyleSheet("""QLineEdit {color: green }""")
        self.addmag1.setStyleSheet("""QDoubleSpinBox {color: red }""")
        self.addmag2.setStyleSheet("""QDoubleSpinBox {color: green }""")
        self.objmag1.setStyleSheet("""QDoubleSpinBox {color: red }""")
        self.objmag2.setStyleSheet("""QDoubleSpinBox {color: green }""")
        self.objna1.setStyleSheet("""QDoubleSpinBox {color: red }""")
        self.objna2.setStyleSheet("""QDoubleSpinBox {color: green }""")

        # define default values for objectives and update ui elements
        objmag = 20
        objna = 0.7
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
        self.mic1 = Microscope(name="Mic1", objmag=objmag, objna=objna, addmag=addmag)

        self.mic2 = Microscope(name="Mic2", objmag=objmag, objna=objna, addmag=addmag)

        self.emwl_value = emwl
        self.emwl.setValue(emwl)
        self.phf1_value = phf
        self.phf1.setValue(phf)
        self.sampling_value = sampling
        self.nyq.setValue(sampling)

        # define camera types
        camera_types = ["CCD", "EM-CCD", "CMOS"]

        # define defaults for camera 1
        name1 = "cam1"
        type1 = "CCD"
        gain1 = 1
        bin1 = 1
        qe1 = 0.74
        pixsize1 = 4.54
        readout1 = 6.5
        dark1 = 0.06
        cic1 = 0.0

        # define defaults for camera 2
        name2 = "cam2"
        type2 = "EM-CCD"
        gain2 = 150
        bin2 = 1
        qe2 = 0.94
        pixsize2 = 12.0
        readout2 = 130
        dark2 = 0.0003
        cic2 = 0.006

        # initialize tow cameras with default values
        self.cam1 = Camera(
            name=name1,
            qe=qe1,
            pixsize=pixsize1,
            binning=bin1,
            cameratype=type1,
            emgain=gain1,
            readout=readout1,
            dark=dark1,
            cic=cic1,
        )

        # update the UI elements with the choosen defaults
        self.name1.setText(name1)
        self.qe1.setValue(qe1)
        self.type1.setCurrentIndex(camera_types.index(type1))
        self.bin1.setCurrentIndex(bin1 - 1)
        self.pixsize1.setValue(pixsize1)
        self.readnoise1.setValue(readout1)
        self.emgain1.setValue(gain1)
        self.dark1.setValue(dark1)
        self.cic1.setValue(cic1)

        self.cam2 = Camera(
            name=name2,
            qe=qe2,
            pixsize=pixsize2,
            binning=bin2,
            cameratype=type2,
            emgain=gain2,
            readout=readout2,
            dark=dark2,
            cic=cic2,
        )

        # check the camera type and enable or disable the gain
        self.checktype()

        self.name2.setText(name2)
        self.qe2.setValue(qe2)
        self.type2.setCurrentIndex(camera_types.index(type2))
        self.bin2.setCurrentIndex(bin2 - 1)
        self.pixsize2.setValue(pixsize2)
        self.readnoise2.setValue(readout2)
        self.emgain2.setValue(gain2)
        self.dark2.setValue(dark2)
        self.cic2.setValue(cic2)

        # adapt the noise factor and readout noise
        self.cam1 = adapt_noise_readout(self.cam1)
        self.cam2 = adapt_noise_readout(self.cam2)
        self.noisef1.setText(str(self.cam1.nf))
        self.noisef2.setText(str(self.cam2.nf))

        # calculate the values for both cameras
        self.cp1, self.cp2 = calc_values(
            self.cam1,
            self.cam2,
            self.mic1,
            self.mic2,
            emwl=self.emwl_value,
            phf=self.phf1_value,
            sampling=self.sampling_value,
        )

        print("Cameras initialized and values calculated.")

        # update ui
        self.phf2.setText(str(self.cp2["flux"]))
        self.sizef1.setText("1.00")
        self.sizef2.setText(str(self.cp2["corrf_pixarea"]))

        # update values for the pixel sizes
        self.piximage1.setText(str(self.cp1["piximage"]))
        self.piximage2.setText(str(self.cp2["piximage"]))
        self.pixrequired1.setText(str(self.cp1["req_pixsize"]))
        self.pixrequired2.setText(str(self.cp2["req_pixsize"]))

        # configure plot
        self.MplWidget.canvas.axes.set_title("Camera SNR Plot", size=18, weight="bold")
        self.MplWidget.canvas.axes.set_xlabel(
            "Photons / Pixel / Frame", size=14, weight="bold"
        )
        self.MplWidget.canvas.axes.set_ylabel("SNR Ratio", size=14, weight="bold")
        self.MplWidget.canvas.axes.grid(True, linestyle="--")
        self.MplWidget.canvas.axes.set_xlim(0, 200)
        self.MplWidget.canvas.axes.set_ylim(0, 10)

        # plot SNR curves (returns a tuple of line objects, thus the comma)
        (self.snr1_curve,) = self.MplWidget.canvas.axes.plot(
            self.cp1["phf"], self.cp1["snr"], "r-", lw=4, label="SNR 1"
        )
        (self.snr2_curve,) = self.MplWidget.canvas.axes.plot(
            self.cp1["phf"], self.cp2["snr"], "g-", lw=4, label="SNR 2"
        )

        # plot indicator lines (returns a tuple of line objects, thus the comma)
        (self.indicator1_line,) = self.MplWidget.canvas.axes.plot(
            self.cp1["phindx"], self.cp1["phindy"], "r--", lw=3, label="PH 1"
        )
        (self.indicator2_line,) = self.MplWidget.canvas.axes.plot(
            self.cp2["phindx"], self.cp2["phindy"], "g--", lw=3, label="PH 2"
        )
        self.MplWidget.canvas.axes.legend()

        self.update_plot()

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

        # connect photon flux  value
        self.phf1.valueChanged.connect(self.change_flux)

    # check camera type and selected gain
    def checktype(self) -> None:
        """
        Checks the type of the camera and adjusts the settings accordingly.

        This method checks the type of both cameras (cam1 and cam2). If the camera type is CMOS or CCD, it disables the EM gain,
        sets the EM gain value to 1, and sets the CIC (Clock Induced Charge) to zero. If the camera type is EM-CCD, it enables the EM gain.

        """

        # disable EM gain if CMOS or CCD
        if self.cam1.cameratype == "CMOS" or self.cam1.cameratype == "CCD":
            # disable the spinbox
            self.emgain1.setDisabled(True)
            # set emgain value for cam1 = 1
            self.cam1.emgain = 1
            # set the value for the spinbox = 1
            self.emgain1.setValue(1)
            # set CIC to zero
            self.cic1.setValue(0.0)
            self.cam1.cic = 0.0
        elif self.cam1.cameratype == "EM-CCD":
            # enable the spinbox
            self.emgain1.setEnabled(True)

        if self.cam2.cameratype == "CMOS" or self.cam2.cameratype == "CCD":
            self.emgain2.setDisabled(True)
            self.cam1.emgain = 1
            self.emgain2.setValue(1)
            self.cic2.setValue(0.0)
            self.cam2.cic = 0.0
        elif self.cam2.cameratype == "EM-CCD":
            self.emgain2.setEnabled(True)

    # modify plot

    def change_scale(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the scale of the plot in the MainWindow.

        This method changes the range for both the x and y axes of the plot according to the values provided by the user
        through the xscale_min, xscale_max, yscale_min, and yscale_max input fields. After changing the scale, it updates the plot.

        """
        # change the range for both axis
        self.MplWidget.canvas.axes.set_xlim(
            self.xscale_min.value(), self.xscale_max.value()
        )
        self.MplWidget.canvas.axes.set_ylim(
            self.yscale_min.value(), self.yscale_max.value()
        )

        # update the plot
        self.MplWidget.canvas.draw()

    # change camera parameters

    def change_binning(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the binning values of the cameras and updates the plot.

        This method changes the binning values of both cameras (cam1 and cam2) according to the current index of the binning
        combo boxes (bin1 and bin2). It then adapts the noise factor and readout noise for both cameras, updates the noise factor
        text fields (noisef1 and noisef2), and updates the plot.

        """

        # change binning values
        self.cam1.binning = self.bin1.currentIndex() + 1
        self.cam2.binning = self.bin2.currentIndex() + 1

        # adapt the noise factor and readout noise
        self.cam1 = adapt_noise_readout(self.cam1)
        self.cam2 = adapt_noise_readout(self.cam2)
        self.noisef1.setText(str(self.cam1.nf))
        self.noisef2.setText(str(self.cam2.nf))

        # update the plot and redraw
        self.update_plot()

    def change_qe(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the quantum efficiency (qe) values of the cameras and updates the plot.

        This method changes the quantum efficiency values of both cameras (cam1 and cam2) according to the current values
        of the quantum efficiency input fields (qe1 and qe2). It then updates the plot.

        """
        # change the qe values values
        self.cam1.qe = self.qe1.value()
        self.cam2.qe = self.qe2.value()

        # update the plot and redraw
        self.update_plot()

    def change_pix(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the pixel size values of the cameras and updates the plot.

        This method changes the pixel size values of both cameras (cam1 and cam2) according to the current values
        of the pixel size input fields (pixsize1 and pixsize2). It then updates the plot.

        """

        # change the pixel size values
        self.cam1.pixsize = self.pixsize1.value()
        self.cam2.pixsize = self.pixsize2.value()

        # update the plot and redraw
        self.update_plot()

    def change_readoutnoise(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the readout noise values of the cameras and updates the plot.

        This method changes the readout noise values of both cameras (cam1 and cam2) according to the current values
        of the readout noise input fields (readnoise1 and readnoise2). It then adapts the noise factor and readout noise for both cameras,
        updates the noise factor text fields (noisef1 and noisef2), and updates the plot.

        """
        # change the readout noise
        self.cam1.readout = self.readnoise1.value()
        self.cam2.readout = self.readnoise2.value()

        # adapt the noise factor and readout noise
        self.cam1 = adapt_noise_readout(self.cam1)
        self.cam2 = adapt_noise_readout(self.cam2)
        self.noisef1.setText(str(self.cam1.nf))
        self.noisef2.setText(str(self.cam2.nf))

        # update the plot and redraw
        self.update_plot()

    def change_dark(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the dark current values of the cameras and updates the plot.

        This method changes the dark current values of both cameras (cam1 and cam2) according to the current values
        of the dark current input fields (dark1 and dark2). It then updates the plot.

        """
        # change the readout noise
        self.cam1.dark = self.dark1.value()
        self.cam2.dark = self.dark2.value()

        # update the plot and redraw
        self.update_plot()

    def change_cic(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the clock-induced charge (CIC) values of the cameras and updates the plot.

        This method changes the CIC values of both cameras (cam1 and cam2) according to the current values
        of the CIC input fields (cic1 and cic2). It then updates the plot.

        """
        # change the readout noise
        self.cam1.cic = self.cic1.value()
        self.cam2.cic = self.cic2.value()

        # update the plot and redraw
        self.update_plot()

    def change_gain(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the gain values of the cameras and updates the plot.

        This method changes the gain values of both cameras (cam1 and cam2) according to the current values
        of the gain input fields (emgain1 and emgain2). It then updates the plot.

        """
        # change the camera gain
        self.cam1.emgain = self.emgain1.value()
        self.cam2.emgain = self.emgain2.value()

        # update the plot and redraw
        self.update_plot()

    def change_type(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the camera type and updates the plot.

        This method changes the camera type of both cameras (cam1 and cam2) according to the current text
        of the camera type combo boxes (type1 and type2). If the camera type is CCD or CMOS, it sets the EM gain to 1.
        It then checks the camera type and adjusts the UI, adapts the noise factor and readout noise for both cameras,
        updates the noise factor text fields (noisef1 and noisef2), and updates the plot.

        """
        # change the camera type
        self.cam1.cameratype = self.type1.currentText()
        self.cam2.cameratype = self.type2.currentText()

        if self.cam1.cameratype == "CCD" or self.cam1.cameratype == "CMOS":
            self.emgain1.setValue(1)
        if self.cam2.cameratype == "CCD" or self.cam2.cameratype == "CMOS":
            self.emgain2.setValue(1)

        # the the camera type and adjust UI
        self.checktype()

        # adapt the noise factor and readout noise
        self.cam1 = adapt_noise_readout(self.cam1)
        self.cam2 = adapt_noise_readout(self.cam2)
        self.noisef1.setText(str(self.cam1.nf))
        self.noisef2.setText(str(self.cam2.nf))

        # update the plot and redraw
        self.update_plot()

    # change optics

    def change_addmag(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the additional magnification values of the cameras and updates the plot.

        This method changes the additional magnification values of both cameras (mic1 and mic2) according to the current values
        of the additional magnification input fields (addmag1 and addmag2). It then updates the plot.

        """
        # change the readout noise
        self.mic1.addmag = self.addmag1.value()
        self.mic2.addmag = self.addmag2.value()

        # update the plot and redraw
        self.update_plot()

    def change_objmag(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the objective magnification values of the cameras and updates the plot.

        This method changes the objective magnification values of both cameras (mic1 and mic2) according to the current values
        of the objective magnification input fields (objmag1 and objmag2). It then updates the plot.

        """
        # change the readout noise
        self.mic1.objmag = self.objmag1.value()
        self.mic2.objmag = self.objmag2.value()

        # update the plot and redraw
        self.update_plot()

    def change_objna(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the objective numerical aperture values of the cameras and updates the plot.

        This method changes the objective numerical aperture values of both cameras (mic1 and mic2) according to the current values
        of the objective numerical aperture input fields (objna1 and objna2). It then updates the plot.

        """
        # change the readout noise
        self.mic1.objna = self.objna1.value()
        self.mic2.objna = self.objna2.value()

        # update the plot and redraw
        self.update_plot()

    def change_sampling(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the sampling value and updates the plot.

        This method changes the sampling value according to the current value of the Nyquist sampling input field (nyq).
        It then updates the plot.

        """
        # change the sampling value
        self.sampling_value = self.nyq.value()

        # update the plot and redraw
        self.update_plot()

    def change_emwl(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the emission wavelength value and updates the plot.

        This method changes the emission wavelength value according to the current value of the emission wavelength input field (emwl).
        It then updates the plot.

        """
        # change the readout noise
        self.emwl_value = self.emwl.value()

        # update the plot and redraw
        self.update_plot()

    # change the photon flux

    def change_flux(self: QtWidgets.QMainWindow) -> None:
        """
        Changes the photon flux value and updates the plot.

        This method changes the photon flux value according to the current value of the photon flux input field (phf1).
        It then updates the plot.

        """
        # change the readout noise
        self.phf1_value = self.phf1.value()

        # update the plot and redraw
        self.update_plot()

    # update the plot and UI

    def update_plot(self):
        """
        Updates the plot and UI.

        This method recalculates the values, updates the line data for the SNR curves and the indicator lines, updates the size factors
        and photon flux, updates pixel sizes, sets the background color of the pixel size in the UI, and updates the whole plot.

        """
        # recalculate the values
        self.cp1, self.cp2 = calc_values(
            self.cam1,
            self.cam2,
            self.mic1,
            self.mic2,
            emwl=self.emwl_value,
            phf=self.phf1_value,
            sampling=self.sampling_value,
        )

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

        # update the size factors and photon flux
        self.sizef2.setText(str(self.cp2["corrf_pixarea"]))
        self.phf2.setText(str(self.cp2["flux"]))

        # update pixel sizes
        self.piximage1.setText(str(self.cp1["piximage"]))
        self.piximage2.setText(str(self.cp2["piximage"]))
        self.pixrequired1.setText(str(self.cp1["req_pixsize"]))
        self.pixrequired2.setText(str(self.cp2["req_pixsize"]))

        # setting background color ofr the pixelsize in the UI
        if self.cp1["piximage"] > self.cp1["req_pixsize"]:
            # if the pixel size is to big set it to orange
            self.piximage1.setStyleSheet("""QLineEdit { background-color: orange;}""")
        else:
            # if the pixel size is to big set it to gree
            self.piximage1.setStyleSheet(
                """QLineEdit { background-color: lightgreen;}"""
            )

        if self.cp2["piximage"] > self.cp2["req_pixsize"]:
            self.piximage2.setStyleSheet("""QLineEdit { background-color: orange;}""")
        else:
            self.piximage2.setStyleSheet(
                """QLineEdit { background-color: lightgreen;}"""
            )

        # update the whole plot
        self.MplWidget.canvas.draw()
        self.MplWidget.canvas.flush_events()


class Camera:
    def __init__(
        self,
        name: str = "Camera1",
        qe: float = 0.75,
        pixsize: float = 4.54,
        binning: int = 1,
        cameratype: str = "CCD",
        emgain: int = 1,
        readout: float = 1.0,
        noisefactor: float = 1.0,
        dark: float = 0.005,
        cic: float = 0.0,
    ) -> None:
        """
        Initializes a new instance of the Camera class.

        Args:
            name (str, optional): Name of the camera. Defaults to "Camera1".
            qe (float, optional): Quantum efficiency of the camera. Defaults to 0.75.
            pixsize (float, optional): Physical pixel size of the camera. Defaults to 4.54.
            binning (int, optional): Binning factor. Defaults to 1.
            cameratype (str, optional): Type of camera - CCD, CMOS or EM-CCD. Defaults to "CCD".
            emgain (int, optional): EM-Gain, will be set to 1 for CCD or CMOS. Defaults to 1.
            readout (float, optional): Readout noise. Defaults to 1.0.
            noisefactor (float, optional): Noise factor for EM-CCD. Defaults to 1.0.
            dark (float, optional): Dark current. Defaults to 0.005.
            cic (float, optional): Clock-induced charge. Defaults to 0.0.

        """

        # allowed types are
        if not cameratype in ["CCD", "CMOS", "EM-CCD"]:
            cameratype = "CCD"
            print("Specified CameraType is not valid. Use CCD as fallback")

        # store all the parameters
        self.name = name
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
    def __init__(
        self,
        name: str = "Mic1",
        objmag: float = 20.0,
        objna: float = 0.95,
        addmag: float = 1.0,
    ) -> None:
        """
        Initializes a new instance of the Microscope class.

        Args:
            name (str, optional): Name of the objective, e.g., the respective "arm" of the detection system. Defaults to "Mic1".
            objmag (float, optional): Magnification factor of the microscope. Defaults to 20.0.
            objna (float, optional): Numerical aperture of the objective. Defaults to 0.95.
            addmag (float, optional): Additional magnification in front of the camera (C-mount etc.). Defaults to 1.0.

        """

        self.mic1 = name
        self.objmag = objmag
        self.objna = objna
        self.addmag = addmag


def calc_values(
    cam1: type[Camera],
    cam2: type[Camera],
    mic1: type[Microscope],
    mic2: type[Microscope],
    emwl: int = 520,
    phf: int = 50,
    sampling: float = 2.0,
) -> Tuple[Dict, Dict]:
    """
    Calculates various values related to two cameras and two microscopes.

    Args:
        cam1 (type[Camera]): The first camera.
        cam2 (type[Camera]): The second camera.
        mic1 (type[Microscope]): The first microscope.
        mic2 (type[Microscope]): The second microscope.
        emwl (int, optional): Emission wavelength. Defaults to 520.
        phf (int, optional): Photon flux. Defaults to 50.
        sampling (float, optional): Sampling rate. Defaults to 2.0.

    Returns:
        Tuple[Dict, Dict]: Two dictionaries containing calculated values for each camera.

    """

    cp1 = {}
    cp2 = {}
    cp1["flux"] = phf

    # pixel size in image plane incl. binning
    cp1["piximage"] = float(
        np.round(cam1.pixsize * cam1.binning / (mic1.objmag * mic1.addmag), 3)
    )
    cp2["piximage"] = float(
        np.round(cam2.pixsize * cam2.binning / (mic2.objmag * mic2.addmag), 3)
    )

    # required pixel since in image to fulfil Nyquist
    cp1["req_pixsize"] = float(
        np.round(0.61 * (emwl / 1000) / (sampling * mic1.objna), 3)
    )
    cp2["req_pixsize"] = float(
        np.round(0.61 * (emwl / 1000) / (sampling * mic2.objna), 3)
    )

    # correction factor for pixel area
    cp2["corrf_pixarea"] = 1.00
    cp2["corrf_pixarea"] = float(
        np.round((cp2["piximage"] ** 2) / (cp1["piximage"] ** 2), 2)
    )

    # create ph vector containing the number of detected photons and use for both cameras
    cp1["phf"] = np.arange(0, 500, 1, dtype=np.int16)

    # calculation of SNR including CIC - Clock Induced Charge
    cp1["snr"] = (
        cam1.qe
        * cp1["phf"]
        / np.sqrt(
            cam1.nf**2 * (cam1.qe * cp1["phf"] + cam1.dark**2 + cam1.cic**2)
            + (cam1.readout_mod**2 / cam1.emgain**2)
        )
    ).astype(float)
    cp2["snr"] = (
        cam2.qe
        * cp1["phf"]
        / np.sqrt(
            cam2.nf**2 * (cam2.qe * cp1["phf"] + cam2.dark**2 + cam2.cic**2)
            + (cam2.readout_mod**2 / cam2.emgain**2)
        )
    ).astype(float)

    # calculate values for photon indicators
    cp2["flux"] = (np.round(cp1["flux"] * cp2["corrf_pixarea"], 0)).astype(int)

    #  calculate explicit SNR values
    cp1["snr_value"] = (
        (cam1.qe * cp1["flux"])
        / np.sqrt(
            cam1.nf**2 * (cam1.qe * cp1["flux"] + cam1.dark**2 + cam1.cic**2)
            + (cam1.readout_mod**2 / cam1.emgain**2)
        )
    ).astype(float)
    cp2["snr_value"] = (
        (cam2.qe * cp2["flux"])
        / np.sqrt(
            cam2.nf**2 * (cam2.qe * cp2["flux"] + cam2.dark**2 + cam2.cic**2)
            + (cam2.readout_mod**2 / cam2.emgain**2)
        )
    ).astype(float)

    cp1["phindx"] = np.array([cp1["flux"], cp1["flux"], 0])
    cp1["phindy"] = np.array([0, cp1["snr_value"], cp1["snr_value"]])

    cp2["phindx"] = np.array([cp2["flux"], cp2["flux"], 0])
    cp2["phindy"] = np.array([0, cp2["snr_value"], cp2["snr_value"]])

    return cp1, cp2


def adapt_noise_readout(cam: Camera) -> Camera:
    """
    Adapts the noise readout of a camera based on its type.

    This function modifies the camera's noise factor, gain, clock-induced charge, and readout noise
    depending on whether it's a CCD, EM-CCD, or CMOS camera.

    Args:
        cam (Camera): The camera whose noise readout is to be adapted.

    Returns:
        Camera: The same camera object with its properties modified based on its type.

    """

    # adjust noise factor due to CCD type
    if cam.cameratype == "CCD":
        # reset noise factor and gain in case of an normal CCD
        cam.nf = 1.0
        cam.emgain = 1
        cam.cic = 0.0
        cam.readout_mod = cam.readout

    elif cam.cameratype == "EM-CCD":
        cam.nf = 1.41
        cam.readout_mod = cam.readout

    # adapt the readout noise if camera is an CMOS
    if cam.cameratype == "CMOS":
        cam.readout_mod = cam.readout * np.sqrt(cam.binning)

    return cam


def main():
    app = QtWidgets.QApplication(sys.argv)
    main = MainWindow()
    main.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
