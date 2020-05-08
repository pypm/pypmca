
from sys import version_info

if version_info.major >= 3:
    from pypm.Population import Population
    from pypm.Model import Model
    from pypm.Delay import Delay
    from pypm.Parameter import Parameter
    from pypm.Multiplier import Multiplier
    from pypm.Propagator import Propagator
    from pypm.Splitter import Splitter
    from pypm.Adder import Adder
    from pypm.Subtractor import Subtractor
    from pypm.Chain import Chain
    from pypm.Modifier import Modifier
    from pypm.Injector import Injector
    from pypm.Ensemble import Ensemble
else:
    raise NotImplementedError('pypm requires Python 3 or higher. Please upgrade!')

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
