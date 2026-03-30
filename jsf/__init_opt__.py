import random
import math
from typing import Any, Callable, Dict, List, NewType, Tuple, Union
from jsf.types import Time, SystemState, CompartmentValue, Trajectory
from jsf import exact
from jsf import sbml

def read_sbml(sbml_xml: str):
    """
    Read an SBML file and return the initial state, rates, and
    stoichiometric matrix.

    Args:
        sbml_xml: The SBML file to read.

    Returns:
        x0: Initial state.
        rates: Function that computes reaction rates.
        stoich: Stoichiometry matrix.
    """
    return sbml.read_sbml(sbml_xml)


def jsf(x0: SystemState, rates, stoich, t_max, **kwargs) -> Trajectory:
    """Generates a sample from the JSF process.

    Args:
        x0: The initial state of the system.
        rates: A function that takes the current state and time and
            returns the rates of each reaction.
        stoich: A dictionary containing the stoichiometry of the
            system.
        t_max: The final time of the simulation.
        **kwargs: A dictionary containing the simulation options.

    Returns:
        A list containing the time series of the state of the system.

    Raises:
        RuntimeError: If the requested method is not implemented.
    """
    method = kwargs['method']
    config = kwargs['config']
    if method is None or method == 'exact':
        return exact.JumpSwitchFlowExact(x0, rates, stoich, t_max, config)
    elif method == 'operator-splitting':
        return JumpSwitchFlowSimulator(x0, rates, stoich, t_max, config)
    else:
        raise RuntimeError(f"Requested method is bonkers {method}")


def JumpSwitchFlowSimulator(
        x0: SystemState,
        rates: Callable[[SystemState, Time], List[float]],
        stoich: Dict[str, Any],
        t_max: Time,
        options: Dict[str, Any]) -> Trajectory:
    """
    Simulate a jump-switch-flow process using the operator splitting method.

    Args:
        x0: The initial state of the system.
        rates: A function that takes the current state and time and
            returns the rates of each reaction.
        stoich: A dictionary containing the stoichiometry of the system.
        t_max: The final time of the simulation.
        options: A dictionary containing the simulation options.

    Returns:
        A list containing the time series of the state of the system.
    """
    # NOTE: In the nu-matrix each row is a reaction and each column
    # describes the net change in that species. One rate per reaction.
    nu = stoich["nu"]               # type: List[List[float]]
    nRates = len(nu)
    nCompartments = len(nu[0])
    nuReactant = stoich["nuReactant"]
    dt = options["dt"]              # type: Time

    # Accept 0/1 integers as well as booleans for backwards compatibility
    EnforceDo = [ed != 0 for ed in options["EnforceDo"]]
    SwitchingThreshold = options["SwitchingThreshold"]  # type: List[int]

    # Determine initial discrete/continuous regime for each compartment
    DoDisc = [x <= thresh and x == round(x)
              for x, thresh in zip(x0, SwitchingThreshold)]

    # Build compartment-in-reaction membership masks.
    # compartInNu[r][i] is True if species i participates in reaction r
    # (either through nu or nuReactant).
    switching_type = options.get('SwitchingType')
    if switching_type is None or switching_type == 'default':
        NuComp    = [[v != 0 for v in row] for row in nu]
        ReactComp = [[v != 0 for v in row] for row in nuReactant]
        compartInNu = [[a or b for a, b in zip(rn, rr)]
                       for rn, rr in zip(NuComp, ReactComp)]
    elif switching_type == 'generous':
        compartInNu = [[v != 0 for v in row] for row in nu]

    # frozenReaction[r] = True means reaction r is handled stochastically (jump),
    # not as continuous flow.
    frozenReaction = [False] * nRates
    for idx in range(nCompartments):
        for reactionIdx in range(nRates):
            if DoDisc[idx] and compartInNu[reactionIdx][idx]:
                frozenReaction[reactionIdx] = True

    # Cumulative integral of propensities and uniform random thresholds
    # used to determine when each reaction fires.
    integralOfFiringTimes = [0.0] * nRates
    randTimes = [random.random() for _ in range(nRates)]
    tauArray  = [Time(0.0)] * nRates

    # Initialise solution storage
    X        = [[x0[i]] for i in range(nCompartments)]
    TauArr   = [Time(0.0)]
    EventType = [0]
    iters    = 0

    AbsT  = Time(0.0)
    ContT = Time(0.0)
    Xprev = x0
    Xcurr = x0

    # newlyDiscCompIndex tracks which compartment most recently switched
    # from continuous to discrete (previously called NewDiscCompartment).
    newlyDiscCompIndex = None
    correctInteger = 0

    _exp = math.exp      # local binding avoids repeated global lookup in tight loop
    _log = math.log
    _floor = math.floor

    while ContT < t_max:

        Dtau  = dt
        Xprev = SystemState([x[iters] for x in X])
        Props = rates(Xprev, ContT)

        # Forward-Euler derivative for continuous compartments
        dXdt = ComputedXdt(Props, nu, frozenReaction, nCompartments)

        # Determine whether any compartment switches regime this step
        (Dtau, correctInteger, DoDisc, frozenReaction,
         NewDoDisc, NewfrozenReaction, newlyDiscCompIndex) = UpdateCompartmentRegime(
            dt, Xprev, Dtau, dXdt, Props, nu, SwitchingThreshold,
            DoDisc, EnforceDo, frozenReaction, compartInNu, nCompartments, nRates)

        # Apply Euler step only to continuous compartments
        Xcurr = SystemState([
            CompartmentValue(X[i][iters] + (0 if DoDisc[i] else Dtau * dXdt[i]))
            for i in range(nCompartments)
        ])

        # If a compartment just became discrete, refresh discrete flags
        OriginalDoDisc = DoDisc[:]
        if correctInteger == 1:
            NewDoDisc, NewfrozenReaction = IsDiscrete(
                Xprev, SwitchingThreshold, DoDisc, EnforceDo,
                frozenReaction, compartInNu, nCompartments, nRates)
            frozenReaction = NewfrozenReaction[:]
            DoDisc         = NewDoDisc[:]

        # --- Stochastic (jump) sub-loop ---
        stayWhile       = any(DoDisc)
        AbsT            = ContT
        DtauContStep    = Dtau
        TimePassed      = Time(0.0)
        firstStayWhileLoop = True

        while stayWhile:

            firstStayWhileLoop = False
            if TimePassed > 0:
                Props = rates(Xcurr, AbsT)

            # Accumulate integral of propensities over this sub-step (trapezoid rule)
            integralStep = ComputeIntegralOfFiringTimes(
                Dtau, Props, rates, Xprev, Xcurr, AbsT)
            integralOfFiringTimes = ArrayPlusAB(
                integralOfFiringTimes,
                ArrayMultiplyAB(integralStep, frozenReaction))

            # If a compartment just became discrete, reset its firing integrals
            if correctInteger == 1:
                for ii in range(nCompartments):
                    if newlyDiscCompIndex == ii and not EnforceDo[ii]:
                        for jj in range(nRates):
                            if compartInNu[jj][ii]:
                                frozenReaction[jj]          = True
                                integralOfFiringTimes[jj]   = 0.0
                                randTimes[jj]               = random.random()

            # Check which reactions have fired (exponential inter-arrival test)
            firedReactions = [
                (0 > (rand - (1 - _exp(-integral)))) and disc
                for rand, integral, disc in
                zip(randTimes, integralOfFiringTimes, frozenReaction)
            ]

            if any(firedReactions):
                tauArray = ComputeFiringTimes(
                    firedReactions, integralOfFiringTimes, randTimes,
                    Props, Dtau, nRates, integralStep)

                if num_non_zero(tauArray) > 0:
                    (Xcurr, Xprev, integralOfFiringTimes, integralStep,
                     randTimes, TimePassed, AbsT, DtauMin, pos) = ImplementFiredReaction(
                        tauArray, integralOfFiringTimes, randTimes, Props, rates,
                        integralStep, TimePassed, AbsT, X, iters, nu, dXdt,
                        OriginalDoDisc, frozenReaction, SwitchingThreshold)

                    # Randomly round any discrete compartment that has a
                    # fractional value, preserving the mean behaviour.
                    for ix, x in enumerate(Xcurr):
                        if (x < SwitchingThreshold[ix]) and (abs(x - round(x)) > 1e-10):
                            x_floor = _floor(x)
                            Xcurr[ix] = (x_floor + 1
                                         if random.uniform(0, 1) < (x - x_floor)
                                         else x_floor)
                            stayWhile = False

                    iters += 1
                    for x_list, xcurr in zip(X, Xcurr):
                        x_list.append(xcurr)
                    TauArr.append(AbsT)
                    EventType.append(pos)

                    Dtau -= DtauMin
                else:
                    stayWhile = False
            else:
                stayWhile = False

            if TimePassed >= DtauContStep:
                stayWhile = False

        # --- End of stochastic sub-loop; finalise continuous step ---
        iters  += 1
        ContT   = Time(ContT + DtauContStep)
        TauArr.append(ContT)
        EventType.append(-1)

        for i in range(len(X)):
            X[i].append(CompartmentValue(
                X[i][iters - 1] + (0 if DoDisc[i] else (DtauContStep - TimePassed) * dXdt[i])
            ))

        # If a compartment switched to discrete this step, round its value
        # and reset the associated firing integrals.
        if correctInteger == 1:
            pos = newlyDiscCompIndex
            X[pos][iters] = round(X[pos][iters])
            for jj in range(nRates):
                if compartInNu[jj][pos]:
                    frozenReaction[jj]          = True
                    integralOfFiringTimes[jj]   = 0.0
                    randTimes[jj]               = random.random()
            newlyDiscCompIndex = None

    return Trajectory((X, TauArr, EventType))


def ComputeFiringTimes(firedReactions, integralOfFiringTimes, randTimes,
                       Props, dt, nRates, integralStep):
    """
    Compute the sub-step time at which each fired reaction occurred,
    using the inverse of the exponential CDF.
    """
    _log = math.log
    tauArray = [0.0] * nRates
    for kk in range(nRates):
        if firedReactions[kk]:
            # Remaining integral needed to reach the firing threshold
            Integral_t0_ti = integralStep[kk] - integralOfFiringTimes[kk]
            Integral       = Integral_t0_ti - _log(1.0 - randTimes[kk])
            tauArray[kk]   = Integral / Props[kk]
    return tauArray


def ImplementFiredReaction(tauArray, integralOfFiringTimes, randTimes, Props,
                           rates, integralStep, TimePassed, AbsT, X, iters,
                           nu, dXdt, OriginalDoDisc, frozenReaction, SwitchingThreshold):
    """
    Apply the earliest-firing reaction, update state and integrals,
    and draw a new random threshold for that reaction channel.
    """
    # Replace zeros (unfired) with inf so min() finds the true earliest
    tauArray  = [float('inf') if tau == 0.0 else tau for tau in tauArray]
    DtauMin   = min(tauArray)
    pos       = tauArray.index(DtauMin)

    TimePassed += DtauMin
    AbsT       += DtauMin

    nCompartments = len(X)
    # Apply deterministic flow + discrete jump to current state
    Xcurr = [
        X[i][iters] + nu[pos][i] + (0 if OriginalDoDisc[i] else DtauMin * dXdt[i])
        for i in range(nCompartments)
    ]
    Xprev = [X[i][iters] for i in range(nCompartments)]

    # Update cumulative integrals: subtract old sub-step, add new sub-step
    integralOfFiringTimes = [
        integral - step * disc
        for integral, step, disc in zip(integralOfFiringTimes, integralStep, frozenReaction)
    ]
    integralStep = ComputeIntegralOfFiringTimes(DtauMin, Props, rates, Xprev, Xcurr, AbsT)
    integralOfFiringTimes = [
        integral + step * disc
        for integral, step, disc in zip(integralOfFiringTimes, integralStep, frozenReaction)
    ]

    # Reset the fired reaction's integral and draw a new threshold
    integralOfFiringTimes[pos] = 0.0
    randTimes[pos]             = random.random()

    return (Xcurr, Xprev, integralOfFiringTimes, integralStep,
            randTimes, TimePassed, AbsT, DtauMin, pos)


def ComputeIntegralOfFiringTimes(Dtau, Props, rates, Xprev, Xcurr, AbsT):
    """
    Integrate propensities over [AbsT, AbsT+Dtau] using the trapezoid rule.
    """
    Props_next = rates(Xcurr, AbsT + Dtau)
    return [Dtau * 0.5 * (p + pn) for p, pn in zip(Props, Props_next)]


def ComputedXdt(
        Props: List[float],
        nu: List[List[float]],
        frozenReaction: List[bool],
        nCompartments: int) -> List[float]:
    """
    Compute the ODE right-hand side for continuous compartments by
    summing propensity-weighted stoichiometry over unfrozen reactions.

    Args:
        Props: Propensities of each reaction.
        nu: Stoichiometry matrix (reactions × species).
        frozenReaction: Which reactions are treated stochastically.
        nCompartments: Number of species/compartments.
    """
    nProps = len(Props)
    return [
        sum(Props[i] * nu[i][j]
            for i in range(nProps) if not frozenReaction[i])
        for j in range(nCompartments)
    ]


def UpdateCompartmentRegime(dt, Xprev, Dtau, dXdt, Props, nu, SwitchingThreshold,
                            DoDisc, EnforceDo, frozenReaction, compartInNu,
                            nCompartments, nRates):
    """
    Check whether any continuous compartment will cross the switching
    threshold within this time step. If so, shorten the step so that
    the crossing happens at the step boundary.
    """
    NewDoDisc, NewfrozenReaction = IsDiscrete(
        Xprev, SwitchingThreshold, DoDisc, EnforceDo,
        frozenReaction, compartInNu, nCompartments, nRates)

    correctInteger     = 0
    newlyDiscCompIndex = None

    # Projected Euler increments for continuous compartments only
    x_step = [0 if isDisc else Dtau * dxi
              for dxi, isDisc in zip(dXdt, NewDoDisc)]

    # Detect any continuous compartment about to fall below its threshold
    crossing = any(
        (x + dxi <= thresh) and not isDisc
        for x, isDisc, thresh, dxi in
        zip(Xprev, NewDoDisc, SwitchingThreshold, x_step)
    )

    if crossing:
        possible_Dtau = [dt]
        for i, (x, isDisc, thresh, dxi) in enumerate(
                zip(Xprev, NewDoDisc, SwitchingThreshold, x_step)):
            if (not isDisc) and (x + dxi <= thresh):
                rounded_x = max(math.ceil(x + x_step[i]), thresh)
                possible_Dtau.append(abs((rounded_x - x) / dXdt[i]))

        if len(possible_Dtau) > 1:
            Dtau = min(possible_Dtau)
            pos  = possible_Dtau.index(Dtau) - 1  # -1 to skip the initial dt entry
            newlyDiscCompIndex = pos
            correctInteger     = 1
        else:
            frozenReaction = NewfrozenReaction
            DoDisc         = NewDoDisc
    else:
        frozenReaction = NewfrozenReaction
        DoDisc         = NewDoDisc

    return (Dtau, correctInteger, DoDisc, frozenReaction,
            NewDoDisc, NewfrozenReaction, newlyDiscCompIndex)


def IsDiscrete(X, SwitchingThreshold, DoDisc, EnforceDo, frozenReaction,
               compartInNu, nCompartments, nRates):
    """
    Determine which compartments should be treated as discrete based on
    whether their current value is at or below the switching threshold.
    EnforceDo overrides the threshold check for pinned compartments.
    """
    DoDiscTmp = [x <= thresh for x, thresh in zip(X, SwitchingThreshold)]

    # Pinned compartments retain their previous discrete/continuous status
    for idx, enforce in enumerate(EnforceDo):
        if enforce:
            DoDiscTmp[idx] = DoDisc[idx]

    # Only recompute frozenReaction if the discrete flags have changed
    if all(DoDiscTmp[i] == DoDisc[i] for i in range(nCompartments)):
        return DoDiscTmp, frozenReaction[:]

    frozenReactionTmp = [False] * nRates
    for idx in range(nCompartments):
        disc_flag = DoDiscTmp[idx] if not EnforceDo[idx] else DoDisc[idx]
        if disc_flag:
            for reactionIdx in range(nRates):
                if compartInNu[reactionIdx][idx]:
                    frozenReactionTmp[reactionIdx] = True

    return DoDiscTmp, frozenReactionTmp


# --- Helper functions ---

def ArraySubtractAB(ArrayA, ArrayB):
    """Element-wise subtraction of two lists."""
    return [a - b for a, b in zip(ArrayA, ArrayB)]

def ArrayPlusAB(ArrayA, ArrayB):
    """Element-wise addition of two lists."""
    return [a + b for a, b in zip(ArrayA, ArrayB)]

def MatrixSubtractAB(MatrixA, MatrixB):
    """Element-wise subtraction of two 2-D lists."""
    return [[a - b for a, b in zip(r1, r2)]
            for r1, r2 in zip(MatrixA, MatrixB)]

def MatrixPlusAB(MatrixA, MatrixB):
    """Element-wise addition of two 2-D lists."""
    return [[a + b for a, b in zip(r1, r2)]
            for r1, r2 in zip(MatrixA, MatrixB)]

def num_non_zero(array: Union[List[float], List[int], List[Time]]) -> int:
    """Count the number of non-zero elements in a list."""
    return sum(1 for element in array if element != 0)

def MatrixDOTArray(Matrix, Array):
    """Multiply a 2-D list (matrix) by a 1-D list (vector)."""
    n = len(Array)
    return [sum(row[i] * Array[i] for i in range(n)) for row in Matrix]

def ArrayMultiplyAB(ArrayA, ArrayB):
    """Element-wise multiplication of two lists."""
    return [a * b for a, b in zip(ArrayA, ArrayB)]
