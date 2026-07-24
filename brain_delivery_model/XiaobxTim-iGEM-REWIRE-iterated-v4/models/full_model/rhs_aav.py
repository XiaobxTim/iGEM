import numpy as np

from models.pbpk.lymphatic_absorption import lymphatic_absorption_rhs
from models.pbpk.organ_distribution import organ_distribution_rhs
from models.bbb.bbb_transport import bbb_transport_rhs
from models.intracellular.aav_intracellular import aav_intracellular_rhs
from models.editing.module5 import editing_rhs


def rhs_aav(t, y, config, idx):
    """
    Total RHS for current AAV model:
    Module 1: absorption
    Module 2: systemic distribution
    Module 3: BBB transport
    Module 4: intracellular uptake and expression
    Module 5: competitive editing
    """
    dydt = np.zeros_like(y)

    # Module 1
    absorption_params = config["absorption"][config["route"]]
    dydt += lymphatic_absorption_rhs(t, y, absorption_params, idx)

    # Module 2
    distribution_params = config["distribution"]
    dydt += organ_distribution_rhs(t, y, distribution_params, idx)

    # Module 3
    bbb_params = config["bbb"]
    dydt += bbb_transport_rhs(t, y, bbb_params, idx)

    # Module 4
    intracellular_params = config["intracellular"]
    dydt += aav_intracellular_rhs(t, y, intracellular_params, idx)

    # Module 5
    editing_params = config["editing"]
    dydt += editing_rhs(t, y, editing_params, idx)

    return dydt
