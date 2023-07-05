"""
This code provides a Bank class to store items and to include methods like
 inserting a proposal to itself or if a proposal has been covered in this bank.
"""
import bisect
from operator import attrgetter

import numpy as np
import matplotlib.pyplot as plt


class LazyNHoods(object):
    __slots__ = ("seq", "nhood_param")

    def __init__(self, seq, nhood_param):
        self.seq = seq
        self.nhood_param = nhood_param

    def __getitem__(self, idx):
        return getattr(self.seq[idx], self.nhood_param)

    def __len__(self):
        return len(self.seq)


class Bank(object):
    __slots__ = ('nhood_size', 'nhood_param', '_templates', '_nmatch', '_NHoods', '_ax')

    def __init__(self, nhood_size=1.0, nhood_param="x1", if_plot=False):
        self.nhood_size = nhood_size
        self.nhood_param = nhood_param

        self._templates = []
        self._nmatch = 0
        self._NHoods = LazyNHoods(self._templates, self.nhood_param)
        if if_plot:
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.set_xlabel('x', fontsize=24)
            ax.set_ylabel('y', fontsize=24)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)  # TODO: improve these hard-coded settings
            ax.tick_params(labelsize=16)
            ax.grid(which='both', zorder=1, alpha=0.5)
            ax.set_aspect('equal')
            fig.tight_layout()
            self._ax = ax
        else:
            self._ax = None

    def __len__(self):
        return len(self._templates)

    def __iter__(self):
        return iter(self._templates)

    def __repr__(self):
        return repr(self._templates)

    def insort(self, new, prefix):
        ind = bisect.bisect_left(self._NHoods, getattr(new, self.nhood_param))
        self._templates.insert(ind, new)
        if self._ax:
            new.ellipse.set(color='C2')
            self._ax.add_patch(new.ellipse)
            dot_tem = self._ax.plot(*new.params, markersize=1, color='C2', marker='.', zorder=3)
            plt.savefig(f'{prefix}_{self._nmatch:0>6d}_{len(self):0>3d}_0.png')
            new.ellipse.set(color='C0')
            for i in dot_tem:
                i.set(color='C0')
            plt.savefig(f'{prefix}_{self._nmatch:0>6d}_{len(self):0>3d}_1.png')

    def add_from_array(self, arr, tmplt_class):
        newtmplts = [tmplt_class(*i) for i in arr]
        # Mark all templates as seed points
        for template in newtmplts:
            template.is_seed_point = True
        self._templates.extend(newtmplts)
        self._templates.sort(key=attrgetter(self.nhood_param))

    @classmethod
    def from_array(cls, arr, tmplt_class, *args, **kwargs):
        bank = cls(*args, **kwargs)
        bank.add_from_array(arr, tmplt_class)
        return bank

    def covers(self, proposal, max_distance, prefix):
        """
        Return (min_distance, template) where min_distance is either
        (i) the best found distance if min_distance < max_distance or
        (ii) the distance of the first template found with distance >= max_distance.
        template is the Template() object which yields min_distance.
        """
        min_distance = np.inf
        template = None

        # find templates in the bank "near" this tmplt
        prop_nhd = getattr(proposal, self.nhood_param)
        low, high = _find_neighborhood(self._NHoods, prop_nhd, self.nhood_size)
        tmpbank = self._templates[low:high]
        if self._ax:
            ellipse_new = proposal.get_ellipse(max_distance, c='C3')
            dot_new = self._ax.plot(*proposal.params, markersize=1, color='C3', marker='.', zorder=3)
            self._ax.add_patch(ellipse_new)
            if self.nhood_param == 'x1':
                p = self._ax.axvspan(prop_nhd-self.nhood_size, prop_nhd+self.nhood_size,
                                     lw=0, color='C1', alpha=0.2, zorder=2)
            elif self.nhood_param == 'x2':
                p = self._ax.axhspan(prop_nhd-self.nhood_size, prop_nhd+self.nhood_size,
                                     lw=0, color='C1', alpha=0.2, zorder=2)
            else:
                raise NotImplementedError
            plt.savefig(f'{prefix}_{self._nmatch:0>6d}_{len(self):0>3d}.png')

        if tmpbank:
            # sort the bank by its nearness to tmplt
            tmpbank.sort(key=lambda b: abs(getattr(b, self.nhood_param) - prop_nhd))

            line_new = []
            # find and test distancees
            for tmplt in tmpbank:

                self._nmatch += 1
                distance = tmplt.proper_distance(proposal)
                if self._ax:
                    line_new += self._ax.plot([proposal.x1, tmplt.x1], [proposal.x2, tmplt.x2],
                                              color='C4', lw=1, zorder=3)
                    # line_new = self._ax.plot([proposal.x1, tmplt.x1], [proposal.x2, tmplt.x2],
                    #                          color='C4', lw=1, zorder=3)
                    # plt.savefig(f'{prefix}_{self._nmatch:0>7d}_{len(self):0>3d}.png')
                    # line_new[0].remove()

                # record distance and template params for largest distance
                if distance < min_distance:
                    min_distance = distance
                    template = repr(tmplt)

                if distance < max_distance:
                    break
            if self._ax:
                plt.savefig(f'{prefix}_{self._nmatch:0>6d}_{len(self):0>3d}.png')
                for i in line_new:
                    i.remove()

        if self._ax:
            ellipse_new.remove()
            dot_new[0].remove()
            p.remove()

        return min_distance, template


def _find_neighborhood(tmplt_locs, prop_loc, nhood_size=0.25):
    """
    Return the min and max indices of templates that cover the given
    template at prop_loc within a parameter difference of nhood_size.
    tmplt_locs should be a sequence of neighborhood values in sorted order.
    """
    low_ind = bisect.bisect_left(tmplt_locs, prop_loc - nhood_size)
    high_ind = bisect.bisect_right(tmplt_locs, prop_loc + nhood_size)
    return low_ind, high_ind
