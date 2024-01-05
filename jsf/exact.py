import random
import math
from typing import Any, Callable, Dict, List, NewType, Tuple, Union
from jsf.types import Time, SystemState, CompartmentValue, Trajectory

JumpClock = NewType('JumpClock', float)
ExtendedState = NewType('ExtendedState',
                        Tuple[SystemState, List[JumpClock], Time])


def JumpSwitchFlowExact(
        x0: SystemState,
        rates: Callable[[SystemState, Time], List[float]],
        stoich: Dict[str, Any],
        t_max: Time,
        options: Dict[str, Any]) -> Trajectory:
    return Trajectory(([], []))


def _update(
        x0: ExtendedState,
        delta_time: Time,
        rates: Callable[[SystemState, Time], List[float]],
        stoich: Dict[str, Any],
        options: Dict[str, Any]) -> ExtendedState:
    return ExtendedState((SystemState([]), [JumpClock(0.0)], Time(0.0)))


def _jump(
        x0: ExtendedState,
        reaction: int,
        stoich: Dict[str, Any],
        options: Dict[str, Any]) -> ExtendedState:
    return ExtendedState((SystemState([]), [JumpClock(0.0)], Time(0.0)))
