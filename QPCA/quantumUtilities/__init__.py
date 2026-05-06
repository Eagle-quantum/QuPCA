import sys

from . import phase_estimation_builder as _phase_estimation_builder
from . import qram_circuit_builder as _qram_circuit_builder
from . import state_vector_tomography as _state_vector_tomography
from .quantum_utilities import thetas_computation,from_binary_tree_to_qcircuit

StateVectorTomography = _state_vector_tomography.StateVectorTomography
PeCircuitBuilder = _phase_estimation_builder.PeCircuitBuilder
QramBuilder = _qram_circuit_builder.QramBuilder

sys.modules[__name__ + ".Tomography"] = _state_vector_tomography
sys.modules[__name__ + ".tomography"] = _state_vector_tomography
sys.modules[__name__ + ".qPe_Builder"] = _phase_estimation_builder
sys.modules[__name__ + ".qpe_builder"] = _phase_estimation_builder
sys.modules[__name__ + ".qRam_Builder"] = _qram_circuit_builder
sys.modules[__name__ + ".qram_builder"] = _qram_circuit_builder

__all__ = [
	'thetas_computation',
	'from_binary_tree_to_qcircuit',
	'StateVectorTomography',
	'PeCircuitBuilder',
	'QramBuilder',
]