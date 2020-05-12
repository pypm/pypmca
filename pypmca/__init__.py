
from sys import version_info

if version_info.major >= 3:
    from pypmca.Population import Population
    from pypmca.Model import Model
    from pypmca.Delay import Delay
    from pypmca.Parameter import Parameter
    from pypmca.Multiplier import Multiplier
    from pypmca.Propagator import Propagator
    from pypmca.Splitter import Splitter
    from pypmca.Adder import Adder
    from pypmca.Subtractor import Subtractor
    from pypmca.Chain import Chain
    from pypmca.Modifier import Modifier
    from pypmca.Injector import Injector
    from pypmca.Ensemble import Ensemble
else:
    raise NotImplementedError('pypmca requires Python 3 or higher. Please upgrade!')

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
