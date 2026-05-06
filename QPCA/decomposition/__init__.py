import sys

from . import qpca_model as _qpca_model

QPCA = _qpca_model.QPCA

sys.modules[__name__ + ".qpca"] = _qpca_model
sys.modules[__name__ + ".Qpca"] = _qpca_model

__all__ = ["QPCA"]
