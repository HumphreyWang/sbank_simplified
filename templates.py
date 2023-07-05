"""
This code provides Template classes, and generators used to propose trial points in some distributions.
"""
import numpy as np
from numpy.random.mtrand import uniform
from matplotlib.patches import Ellipse


def uniform_points_generator(**constraints):
    """Uniformly generate points in 2d."""
    x1_min, x1_max = constraints.pop('x1')
    x2_min, x2_max = constraints.pop('x2')

    while 1:   # This is inexplicably much faster than "while True"
        x1 = uniform(x1_min, x1_max)
        x2 = uniform(x2_min, x2_max)
        yield x1, x2


def cartesian_uniform_generator(tmplt_class, **constraints):
    """Uniformly generate templates in 2d, here (x1, x2) denotes for (x, y)."""
    for x1, x2 in uniform_points_generator(**constraints):
        yield tmplt_class(x1, x2)


def polar_uniform_generator(tmplt_class, **constraints):
    """Uniformly generate templates in 2d, here (x1, x2) denotes for (r, theta), where r>=0, theta in [0, 2*pi]."""
    for x1, x2 in uniform_points_generator(**constraints):
        yield tmplt_class(x1*np.cos(x2), x1*np.sin(x2))


class BasicTemplate(object):
    """Basic class, i.e. the class corresponds to Cartesian coordinates."""
    __slots__ = ("x1", "x2", "norm", "is_seed_point", "ellipse")
    param_names = ("x1", "x2")
    param_formats = ("%.4f", "%.4f")
    metric = np.array([[1, 0], [0, 1]])
    vals, vecs = np.linalg.eig(np.linalg.inv(metric))
    ang = float(np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0])))

    def __init__(self, x1, x2):
        self.x1 = float(x1)
        self.x2 = float(x2)
        # TODO: x1, x2 attributes may denote different quantities with those in sbank.py (rename them?)
        norm_sq = float(np.linalg.multi_dot([[x1, x2], self.metric, [[x1], [x2]]]))
        self.norm = norm_sq**0.5
        self.ellipse = None

    @property
    def params(self):
        return tuple(getattr(self, k) for k in self.param_names)

    def __repr__(self):
        return "(%s)" % ", ".join(self.param_formats) % self.params

    def proper_distance(self, other):
        dx1, dx2 = other.x1-self.x1, other.x2-self.x2
        dis_sq = float(np.linalg.multi_dot([[dx1, dx2], self.metric, [[dx1], [dx2]]]))
        return dis_sq**0.5

    def get_ellipse(self, max_distance, c):
        """Provide a ellipse patch for drawing figures."""
        if not self.ellipse:
            w, h = 2*max_distance*np.sqrt(self.vals)
            self.ellipse = Ellipse((self.x1, self.x2), width=w, height=h, angle=self.ang, color=c, alpha=0.2)
        return self.ellipse


class ScaledEuclidTemplate(BasicTemplate):
    """The class for a scaled Euclidean coordinates, i.e. we stretch the x1-direction with a factor of 2.
     Actually you can arbitrarily change the metric as long as it is not singular,
     non-diagonal matrix will give a rotate ellipse, e.g. you can override the metric using [[1/4, 1/4], [1/4, 1]].
    """
    metric = np.array([[1/4, 0], [0, 1]])
    vals, vecs = np.linalg.eig(np.linalg.inv(metric))
    ang = float(np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0])))


proposals = {'Cartesian': cartesian_uniform_generator,
             'Polar': polar_uniform_generator,
             'ScaledEuclidean': cartesian_uniform_generator,
             }
coord_frames = {'Cartesian': BasicTemplate,
                'Polar': BasicTemplate,
                'ScaledEuclidean': ScaledEuclidTemplate,
                }
