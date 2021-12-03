# -*- coding: utf-8 -*-

#################################################################
# File        : mplwidget.py
# Version     : 0.0.1
# Author      : sebi06
# Date        : 21.10.2021
#
# Disclaimer: The code is purely experimental. Feel free to
# use it at your own risk.
#
#################################################################

from PyQt5.QtWidgets import *

from matplotlib.backends.backend_qt5agg import FigureCanvas

from matplotlib.figure import Figure


class MplWidget(QWidget):

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.canvas = FigureCanvas(Figure())

        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)

        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.setLayout(vertical_layout)