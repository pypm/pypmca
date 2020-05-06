
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
else:
    raise NotImplementedError('pypm requires Python 3 or higher. Please upgrade!')