import numpy as np

from .bogoliubov import BogoliubovState


class BCSState(BogoliubovState):
    """This class implements a BCS mean-field state.

    BCS mean-field states are characterized by occupancies ($\rho$) being a
    diagonal matrix and thus _isa_ `BogoliubovState`.
    """

    def __init__(self, hfpsi, _is_base_class=False) -> None:
        """
        Args:
            hfpsi: HFPsi object, containing the single-particle wave functions.
            _is_base_class: True if called by super().__init__(), False otherwise.
        """
        super().__init__(hfpsi, _is_base_class=True)
        if not _is_base_class:
            self.rho = np.empty(self.n_total_wf, dtype=np.float64)

    def validate(self):
        """Verify that all conditions for representing a BCS mean-field State
        are satisfied.

        Returns:
            None
        Raises:
            AssertionError: If any of the conditions are not satisfied.
        """
        assert self.rho.shape == (self.n_total_wf,)
        assert np.all(self.rho >= 0)
        assert np.all(self.rho <= 1)
