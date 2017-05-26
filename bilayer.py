import os
import numpy as np
class bilayer():
    """ bilayer object to store data

    Convention for offsets
    ---------------------
    OH, FFA, etc etc """
    def __init__(self, name = "" , apl = 0, apl_std = 0, apt = 0, apt_std =0, height = 0,
            height_std = 0 , offsets = 0 , offsets_std = 0, tilt_angle =0,
            tilt_angle_std = 0, idig = 0, idig_std = 0):
        self._name = name
        self._apl = apl
        self._apl_std = apl_std
        self._apt = apt
        self._apt_std = apt_std
        self._height = height
        self._height_std = height_std
        self._offsets = offsets
        sefl._offsets_std = offsets_std
        self._tilt_angle = tilt_angle
        self._tilt_angle_std = tilt_angle_std
        self._idig = idig
        self._idig_std = idig_std

    @property
    def name(self):
        return self._name

    @property
    def apl(self):
        return self._apl

    @property
    def apl_std(self):
        return self._apl_std

    @property
    def apt(self):
        return self._apt

    @property
    def apt_std(self):
        return self._apt_std

    @property
    def height(self):
        return self._height

    @property
    def height_std(self):
        return self._height_std

    @property
    def offsets(self):
        return self._offsets

    @property
    def offsets_std(self):
        return self._offsets_std


    @property
    def tilt_angle(self):
        return self._tilt_angle

    @property
    def tilt_angle_std(self):
        return self._tilt_angle_std

    @property
    def idig(self):
        return self._idig

    @property
    def idig_std(self):
        return self._idig_std

    @name.setter
    def name(self, name):
        self._name = name

    @apl.setter
    def apl(self, apl):
        self._apl = apl

    @apl_std.setter
    def apl_std(self, apl):
        self._apl_std = apl

    @apt.setter
    def apt(self, apt):
        self._apt = apt

    @apt_std.setter
    def apt_std(self, apt):
        self._apt_std = apt


    @height.setter
    def height(self, height):
        self._height = height

    @height_std.setter
    def height_std(self, height):
        self._height_std = height


    @offsets.setter
    def offsets(self, offsets):
        self._offsets = offsets

    @offsets_std.setter
    def offsets_std(self, offsets):
        self._offsets_std = offsets


    @tilt_angle.setter
    def tilt_angle(self, tilt):
        self._tilt_angle = tilt

    @tilt_angle_std.setter
    def tilt_angle_std(self, tilt):
        self._tilt_angle = tilt

    @idig.setter
    def idig(self,idig):
        self._idig = idig

    @idig_std.setter
    def idig_std(self,idig):
        self._idig_std = idig
