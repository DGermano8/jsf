import random
import math
import copy
from typing import Any, Callable, Dict, List, NewType, Tuple, Union
from dataclasses import dataclass
from jsf.types import Time, SystemState, CompartmentValue, Trajectory


@dataclass
class JumpClock:
    clock: Time
    start_time: Time

    def set_clock(self, new_clock: float) -> None:
        self.clock = Time(new_clock)


ExtendedState = NewType('ExtendedState',
                        Tuple[SystemState, List[bool], List[JumpClock], Time])


def JumpSwitchFlowExact(
        x0: SystemState,
        rates: Callable[[SystemState, Time], List[float]],
        stoich: Dict[str, Any],
        t_max: Time,
        options: Dict[str, Any]) -> Trajectory:
    """
    Simulate a jump-switch flow model numerically by using forward
    Euler. It is exact in the sense that it does not use operator
    splitting.

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
    jump_clocks = [_new_jump_clock(is_j) for is_j in is_jumping]
    ext_state = ExtendedState(
        (x0, is_jumping, jump_clocks, Time(0.0))
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


def _is_jumping(
        x0: SystemState,
        stoich: Dict[str, Any],
        options: Dict[str, Any]) -> List[bool]:
    """
    Determine which reactions are jumping.

    Notes:
        A reaction is jumping if at least one of the reactants or
    products is discrete.
    """
    threshold = options['SwitchingThreshold']
    is_discrete = [x0[i] <= threshold[i] for i in range(len(x0))]
    is_reactant = [[n != 0 for n in r] for r in stoich['nuReactant']]
    is_product = [[n != 0 for n in r] for r in stoich['nuProduct']]
    num_reactions = len(is_reactant)
    num_compartments = len(is_discrete)
    result = [False] * num_reactions

    for rix in range(num_reactions):
        if any([is_discrete[ix] for ix in range(num_compartments)
                if is_reactant[rix][ix] or is_product[rix][ix]]):
            result[rix] = True

    return result


def _new_jump_clock(
        is_jumping: bool) -> JumpClock:
    """
    Generate a new jump clock.

    Notes:
        If the system is not jumping then the jump clock is set to
    infinity. Otherwise, it takes a random value from the uniform
    distribution on [0, 1].
    """
    u = Time(random.uniform(0.0, 1.0))
    return JumpClock(
        clock = u if is_jumping else Time(float('inf')),
        start_time = u if is_jumping else Time(float('inf'))
    )


def _update_state(
        state: SystemState,
        is_jumping: List[bool],
        reaction_rates: List[float],
        nu_matrix,
        delta_time: Time) -> SystemState:
    """
    Update the state of the system with a single forward-Euler step.
    """
    dx_dt = [
        sum(0 if is_jumping[rix] else reaction_rates[rix] * nu_matrix[rix][cix]
            for rix in range(len(reaction_rates)))
        for cix in range(len(state))
    ]
    return SystemState([
        state[cix] + dx_dt[cix] * delta_time
        for cix in range(len(state))
    ])


def _update_jump_clocks(
        curr_time: Time,
        delta_time: Time,
        curr_jump_clocks: List[JumpClock],
        curr_state: SystemState,
        curr_reaction_rates: List[float],
        next_state: SystemState,
        next_reaction_rates: List[float]) -> List[JumpClock]:
    """
    Update the jump clocks using the trapezoidal rule and the current
    and next system state to approximate the reaction rates.

    Notes:
        Because `enumerate` returns a copy of each element in the
    list we need to assign this value back at the end of the loop for
    it to be updated.
    """
    next_jump_clocks = copy.deepcopy(curr_jump_clocks)
    for ix, jc in enumerate(next_jump_clocks):
        integral = 0.5 * delta_time * (next_reaction_rates[ix] + curr_reaction_rates[ix])
        jc.set_clock(jc.clock + (math.exp(-integral) - 1) * (jc.clock + 1 - jc.start_time))
        next_jump_clocks[ix] = jc
    return next_jump_clocks


def _update(
        x0: ExtendedState,
        delta_time: Time,
        rates: Callable[[SystemState, Time], List[float]],
        stoich: Dict[str, Any],
        options: Dict[str, Any]) -> ExtendedState:
    """
    Update the state of the system.

    Args:
        x0: Current state.
        delta_time: Time step.
        rates: Function that computes reaction rates.
        stoich: Stoichiometry matrix.
        options: Simulation options.

    Returns:
        Updated state.
    """
    # TODO Implement this!
    raise RuntimeError('Not implemented')
    return ExtendedState((SystemState([]), [JumpClock(0.0)], Time(0.0)))


def _jump(
        x0: ExtendedState,
        reaction: int,
        stoich: Dict[str, Any],
        options: Dict[str, Any]) -> ExtendedState:
    """
    Jump the system to a new state.

    Args:
        x0: Current state.
        reaction: Reaction to jump to.
        stoich: Stoichiometry matrix.
        options: Simulation options.

    Returns:
        New state.
    """
    # TODO Implement this!
    raise RuntimeError('Not implemented')
    return ExtendedState((SystemState([]), [JumpClock(0.0)], Time(0.0)))
