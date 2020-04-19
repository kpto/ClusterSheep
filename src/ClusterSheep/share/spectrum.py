# -*- coding: utf-8 -*-
"""
Created on 06:40:37 14/11/2017

Author: Paul TO, Ka Po
Contact: kpto@ust.hk
Institude: Hong Kong University of Science and Technology
Department: Division of Biomedical Engineering

Supervisor: Prof. Henry LAM, H. N.
Contact: kehlam@ust.hk
Institude: Hong Kong University of Science and Technology
Department: Department of Chemical and Biomolecular Engineering

Desciption of this module:
"""

# ====BEGIN OF MODULE IMPORT====
import logging
import numpy as np
import matplotlib.pyplot as plt
import pyopenms as om
from pprint import pformat

from ClusterSheep.envr.session import get_session
from ClusterSheep.share.internal_index import Entry
import ClusterSheep.reader.mzxml as mzxml
import ClusterSheep.reader.mzml as mzml
# ====END OF MODULE IMPORT====


# ====BEGIN OF CONSTANT DEFINITION====
# ====END OF CONSTANT DEFINITION====


# ====BEGIN OF GLOBAL VARIABLE DECLARATION====
try:
    session = get_session()
except ImportError:
    err_msg = '\nThis module requires a valid session.'
    logging.error(err_msg)
    raise
# ====END OF GLOBAL VARIABLE DECLARATION====


# ====BEGIN OF CLASS DEFINITION====
class Spectrum(Entry):

    def __init__(self, obj=None):
        self.mz = None
        self.intensity = None
        self.ranked = False
        self.override_iden = None
        type_ = type(obj)
        if type_ is Entry or type_ is Spectrum:
            for item in vars(obj).items():
                setattr(self, *item)
            if type_ is Entry:
                self.mz, self.intensity = self._get_peaks()
        return

    def _get_peaks(self):
        if self.format.lower() == '.mzxml':
            return mzxml.get_peaks(self.get_file_path(), self.offset)
        elif self.format.lower() == '.mzml':
            return mzml.get_peaks(self.get_file_path(), self.offset)

    def is_sorted(self):
        if len(self.mz) < 2:
            return True
        return np.all(self.mz[1:] > self.mz[:-1])

    def sort_by_mz(self, copy=True):
        target = Spectrum(self) if copy else self
        idx = np.argsort(target.mz)
        target.mz = target.mz[idx]
        target.intensity = target.intensity[idx]
        return target

    def clip(self, mz_range=None, copy=True):
        mz_range = mz_range if mz_range else session.config.rt_mz_range.value
        target = Spectrum(self) if copy else self
        if mz_range == (0, 0): return target

        idx = np.logical_and(target.mz >= mz_range[0], target.mz <= mz_range[1])
        target.mz, target.intensity = target.mz[idx], target.intensity[idx]
        return target

    def remove_precursor(self, removal_range=None, true_precursor_mass=None, copy=True):
        removal_range = removal_range if removal_range else session.config.rt_precursor_removal_range.value
        true_precursor_mass = true_precursor_mass if true_precursor_mass else session.config.ii_true_precursor_mass.value
        target = Spectrum(self) if copy else self
        if removal_range == (0.0, 0.0): return target

        precursor_mass = target.precursor_mass / target.precursor_charge if true_precursor_mass else target.precursor_mass
        idx = np.logical_or(target.mz < (removal_range[0] + precursor_mass),
                            target.mz > (removal_range[1] + precursor_mass))
        target.mz, target.intensity = target.mz[idx], target.intensity[idx]
        return target

    def rank_transform(self, num_of_peaks=None, bins_per_th=None, copy=True, normalize=True):
        num_of_peaks = num_of_peaks if num_of_peaks else session.config.rt_num_of_peaks.value
        bins_per_th = bins_per_th if bins_per_th else session.config.rt_bins_per_th.value
        target = Spectrum(self) if copy else self

        if not target.is_sorted(): target.sort_by_mz(copy=False)
        mz = target.mz if target.mz.dtype == np.float32 else target.mz.astype(np.float32)
        intensity = target.intensity if target.intensity.dtype == np.float32 else target.intensity.astype(np.float32)

        if bins_per_th != 1: mz *= bins_per_th
        mz, intensity = target._binning(mz, intensity)

        temp = min(len(mz), num_of_peaks)
        new_intensity = np.arange(1, temp + 1)[::-1]
        idx = np.argsort(intensity)[:-(num_of_peaks + 1):-1]
        intensity[idx] = new_intensity
        idx.sort()

        mz = mz[idx]
        intensity = intensity[idx]
        if normalize:
            intensity /= np.linalg.norm(intensity)

        if temp < num_of_peaks:
            mz = np.append(np.full(num_of_peaks - temp, -1, dtype=np.int32), mz)
            intensity = np.append(np.full(num_of_peaks - temp, 0, dtype=np.float32), intensity)
            
        target.mz, target.intensity = mz, intensity
        target.ranked = True
        return target

    def ranked_dot_product(self, target):
        if self.ranked is not True or target.ranked is not True: return None
        if len(self.mz) != len(target.mz): return None
        num_of_peaks = len(self.mz)
        self_ptr = 0
        target_ptr = 0
        dp = np.float32(0.0)
        while self_ptr < num_of_peaks and target_ptr < num_of_peaks:
            if self.mz[self_ptr] == target.mz[target_ptr]:
                dp += self.intensity[self_ptr] * target.intensity[target_ptr]
                self_ptr += 1
                target_ptr += 1
            elif self.mz[self_ptr] < target.mz[target_ptr]:
                self_ptr += 1
            else:
                target_ptr += 1
        return dp
    
    def verificative_ranked_dp(self, target):
        a = self.clip().remove_precursor(copy=False).rank_transform(copy=False)
        b = target.clip().remove_precursor(copy=False).rank_transform(copy=False)
        return Spectrum.ranked_dot_product(a, b)

    @staticmethod
    def _binning(mz, intensity):
        mz = mz.astype(np.int32)
        new_intensity = np.zeros(mz.max()+1, dtype=np.float32)
        np.add.at(new_intensity, mz, intensity)
        mz = np.unique(mz)
        return mz, new_intensity[mz]

    def _plot(self, ax, color, top=True, verificative=True, show_iden=True, highlight=True):
        target = self if self.is_sorted() else self.sort_by_mz()
        if verificative:
            temp = target.clip().remove_precursor(copy=False).rank_transform(copy=False)
            mz, intensity = temp.mz, temp.intensity
        else:
            mz, intensity = target._binning(target.mz, target.intensity)
        rects = ax.bar(mz, intensity, 1, color=color, edgecolor=color)
        if show_iden or highlight:
            identification = target.get_identification() if not self.override_iden else self.override_iden
            if identification:
                tpp_string = identification.to_tpp_string()
                title = identification.to_string(tpp_string)
                if self.override_iden: title = 'Overridden: ' + title
                if show_iden:
                    y_loc = 0.9 if top else 0.1
                    ax.text(0.5, y_loc, title, horizontalalignment='center',
                            verticalalignment='center', transform=ax.transAxes, fontsize=36)
                if highlight:
                    tspec = om.MSSpectrum()
                    try:
                        aa_sequence = om.AASequence.fromString(tpp_string, False)
                    except Exception:
                        logging.info('String "{}" cannot be understood by pyopenms.AASequence.'.format(tpp_string))
                        aa_sequence = None
                    if aa_sequence:
                        Spectrum.tsg.getSpectrum(tspec, aa_sequence, 1, identification.charge)
                        ions_mz = tspec.get_peaks()[0].astype(np.uint32)
                        annotations = tspec.getStringDataArrays()[0]
                        temp_dict = {}
                        colors = {'b': 'blue',
                                  'y': 'red',
                                  'a': 'green',
                                  'c': 'teal',
                                  'x': 'purple',
                                  'z': 'orange'}
                        for i in range(len(ions_mz)):
                            temp_dict[ions_mz[i]] = annotations[i]
                        for i in range(len(mz)):
                            if mz[i] in temp_dict.keys():
                                ion_name = temp_dict[mz[i]].decode()
                                color = colors[ion_name[0]]
                                rect = rects[i]
                                height = rect.get_height()
                                ax.text(rect.get_x() + rect.get_width() / 2., 1.05 * height, ion_name, fontsize=18,
                                        ha='center', va='center', color=color)
                                rect.set_color(color)
                                rect.set_edgecolor(color)
        return

    def plot(self, against=None, x_limit=(0, 2000), verificative=True, show_iden=True, highlight=True):
        fig, ax = plt.subplots(2, 1, figsize=(23, 11.5), sharex=True)
        color_top = 'silver'
        color_bottom = 'silver'

        ax[0].set_xlim(x_limit)
        self._plot(ax[0], color_top, top=True, verificative=verificative, show_iden=show_iden, highlight=highlight)
        if against:
            against._plot(ax[1], color_bottom, top=False, verificative=verificative, show_iden=show_iden, highlight=highlight)
            ax[1].set_xlim(x_limit)
        ax[1].invert_yaxis()
        fig.subplots_adjust(hspace=0)
        plt.show()
        return

    def get_spectrum(self):
        return Spectrum(session.internal_index[self.internal_id])

    def copy(self):
        return Spectrum(self)

    def __repr__(self):
        return pformat(vars(self))
# ====END OF CLASS DEFINITION====


# ====BEGIN OF CODE====
tsg = om.TheoreticalSpectrumGenerator()
param = tsg.getParameters()
param.setValue('add_metainfo', b'true', '', [])
param.setValue('add_a_ions', b'true', '', [])
param.setValue('add_c_ions', b'true', '', [])
param.setValue('add_x_ions', b'true', '', [])
param.setValue('add_z_ions', b'true', '', [])
tsg.setParameters(param)
Spectrum.tsg = tsg
# ====END OF CODE====


# ====CODE FOR MODULE TEST====
if __name__ == '__main__':
    pass
# ====CODE FOR MODULE TEST====
