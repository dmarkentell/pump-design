"""Motor selection and the shaft diameter at the impeller, from the duty point alone.

These come first because the impeller hub has to fit around the shaft, and the hub sets the
eye flow area that the hydraulic design depends on. They use only the duty and the
statistical efficiency, so there is no circular dependency on the detailed hydraulic model,
and the JavaScript port can apply the identical rule.

Motor: API 610 asks that the driver nameplate rating be at least a percentage of the pump's
rated power, larger for small pumps: 125% up to 22 kW, 115% from 22 to 55 kW, 110% above 55 kW.
The motor is the next standard NEMA rating at or above that.

Shaft at the impeller: sized for the motor's full nameplate torque at a design shear stress of
40 MPa, which leaves room for the keyway and fatigue. drive_end.py then checks fatigue,
deflection and critical speed properly; this is only the starting size.
"""
from __future__ import annotations

import math

from .correlations import eta_overall

G = 9.80665
RHO = 998.2
HP_W = 745.7
NEMA_HP = [1, 1.5, 2, 3, 5, 7.5, 10, 15, 20, 25, 30, 40, 50, 60, 75, 100, 125, 150, 200, 250, 300,
           350, 400, 450, 500]
TAU_DESIGN = 40e6   # Pa


def api_driver_factor(p_kw: float) -> float:
    if p_kw <= 22:
        return 1.25
    if p_kw <= 55:
        return 1.15
    return 1.10


def estimated_power_kw(q: float, head: float, n_rpm: float) -> float:
    nq = n_rpm * math.sqrt(q) / head ** 0.75
    return RHO * G * q * head / eta_overall(q, nq) / 1e3


def motor_hp(q: float, head: float, n_rpm: float) -> float:
    p = estimated_power_kw(q, head, n_rpm)
    need_hp = p * api_driver_factor(p) * 1e3 / HP_W
    return next(h for h in NEMA_HP if h >= need_hp - 1e-9)


def motor_torque(q: float, head: float, n_rpm: float) -> float:
    """Nameplate torque, N m."""
    return motor_hp(q, head, n_rpm) * HP_W / (2 * math.pi * n_rpm / 60)


def shaft_at_impeller(q: float, head: float, n_rpm: float) -> float:
    """Shaft diameter at the impeller fit, rounded up to 5 mm. Metres."""
    t = motor_torque(q, head, n_rpm)
    d = (16 * t / (math.pi * TAU_DESIGN)) ** (1 / 3)
    return math.ceil(d * 1e3 / 5 - 1e-9) * 5 / 1e3


def hub_diameter(q: float, head: float, n_rpm: float) -> float:
    """Impeller hub OD at the eye: 1.6 x the shaft, room for the bore, key and nut."""
    return 1.6 * shaft_at_impeller(q, head, n_rpm)
