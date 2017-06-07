from numpy import pi, sin, sinh, cos, arange, subtract, around, exp, insert,arctan
import matplotlib.pyplot as plt
from config import *
import os

class ElectrostaticSuspension(object):
    """A Electrostatic suspension suspension object. It has the
    following properties:
        width: micro-meters, plate width
        length: micro-meters, plate length
        ...
    """
    def __init__(self, name, width, length, eps_r=0.0):
        self.name = name
        self.width = width
        self.length = length
        self.eps_r = eps_r

    def capacitance(self, amount):
        """Return the capacitance of suspension"""
        return self.capacitance

    def calc_cap_vs_gap(self, gap_list):


    def read_cap_vs_gap(self, gap_list):
        """Return the capacitance versus gap data points:
        GAP1, CAP1, GAP2, CAP2 ... GAP20, CAP20"""

        self.cap_vs_gap = 0
        return self.cap_vs_gap

class AnsysSimulation(object):

    def __init__(self, name, path='./'):
        self.name = name
        self.path = path