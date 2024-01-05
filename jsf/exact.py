import random
import math
from typing import Any, Callable, Dict, List, NewType, Tuple, Union
from jsf.types import Time, SystemState, CompartmentValue, Trajectory

JumpClock = NewType('JumpClock', float)

ExtendedState = NewType('ExtendedState',
                        Tuple[SystemState, List[bool], List[JumpClock], Time])


def JumpSwitchFlowExact(
        x0: SystemState,
        rates: Callable[[SystemState, Time], List[float]],
        stoich: Dict[str, Any],
        t_max: Time,
        options: Dict[str, Any]) -> Trajectory:
    """
    Simulate a jump-switch flow model exactly.

    Args:
        x0: Initial state.
        rates: Function that computes reaction rates.
        stoich: Stoichiometry matrix.
        t_max: Maximum time to simulate.
        options: Simulation options.

    Returns:
        Trajectory of the simulation.
    """
    is_jumping = _is_jumping(x0, stoich, options)
    jump_clocks = map(_new_jump_clock, is_jumping)
    ext_state = ExtendedState(
        (curr_state, is_jumping, curr_jump_clocks, Time(0.0))
    )

    time = ext_state[3]
    times = [time]
    state_history = [ext_state[0]]

    while time < t_max:
        ext_state = _update(ext_state, options['dt'], rates, stoich, options)
        time = ext_state[3]
        times.append(time)
        state_history.append(ext_state[0])

    compartment_histories = list(map(list, zip(*state_history)))
    return Trajectory((compartment_histories, times))


def _is_jumpting(
        x0: SystemState,
        stoich: Dict[str, Any],
        options: Dict[str, Any]) -> List[bool]:
    # TODO Implement this!
    raise RuntimeError('Not implemented')
    return []


def _new_jump_clock(
        is_jumping: bool) -> JumpClock:
    raise RuntimeError('Not implemented')
    rand_jump_clock = float('inf') # TODO Implement this!
    return rand_jump_clock if is_jumping else JumpClock(float('inf'))


def _update(
        x0: ExtendedState,
        delta_time: Time,
        rates: Callable[[SystemState, Time], List[float]],
        stoich: Dict[str, Any],
        options: Dict[str, Any]) -> ExtendedState:
    # TODO Implement this!
    raise RuntimeError('Not implemented')
    return ExtendedState((SystemState([]), [JumpClock(0.0)], Time(0.0)))


def _jump(
        x0: ExtendedState,
        reaction: int,
        stoich: Dict[str, Any],
        options: Dict[str, Any]) -> ExtendedState:
    # TODO Implement this!
    raise RuntimeError('Not implemented')
    return ExtendedState((SystemState([]), [JumpClock(0.0)], Time(0.0)))
