import random
import math
from typing import Any, Callable, Dict, List, NewType, Tuple, Union
from jsf.types import Time, SystemState, CompartmentValue, Trajectory


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
    if kwargs['method'] is None or kwargs['method'] == 'exact':
        # throw an error because the method is not implemented
        raise RuntimeError("Exact method not implemented")
    elif kwargs['method'] == 'operator-splitting':
        result = JumpSwitchFlowSimulator(x0, rates, stoich, t_max, kwargs['config'])
    else:
        raise RuntimeError(f"Requested method is bonkers {kwargs['method']}")

    return result


def JumpSwitchFlowSimulator(
        x0: SystemState,
        rates: Callable[[SystemState, Time], List[float]],
        stoich: Dict[str, Any],
        t_max: Time,
        options: Dict[str, Any]) -> Trajectory:
    """
    Simulate a jump-switch-flow process using the operator splitting
    method.

    Args:
        x0: The initial state of the system.
        rates: A function that takes the current state and time and
            returns the rates of each reaction.
        stoich: A dictionary containing the stoichiometry of the
            system.
        t_max: The final time of the simulation.
        options: A dictionary containing the simulation options.

    Returns:
        A list containing the time series of the state of the system.
    """
    # NOTE In the nu-matrix each row is a reaction and each column
    # describes the number items of that species used in the reaction.
    # There is one rate for each reaction and since there is a column
    # for each species we can use the length of the first column to
    # find the number of species in the system.

    nu = stoich["nu"] # type: List[List[float]]
    nRates = len(nu)
    nCompartments = len(nu[0])
    nuReactant = stoich["nuReactant"] # type: List[List[float]]
    dt = options["dt"] # type: Time

    # This is to enable us to use a Boolean list here while still
    # accepting 0/1 as input. In future versions it would be nice to
    # require a Boolean list but that would break backwards
    # compatibility so we will leave it for now.
    EnforceDo = [(not (ed == 0)) for ed in options["EnforceDo"]]

    SwitchingThreshold = options["SwitchingThreshold"] # type: List[int]

    DoDisc = [(x <= threshold and x==round(x)) for x, threshold in zip(x0, SwitchingThreshold)]

    # identify which compartment is in which reaction:
    NuComp = [[value != 0 for value in row] for row in nu]
    ReactComp = [[value != 0 for value in row] for row in nuReactant]
    compartInNu = [[value != 0 for value in row] for row in MatrixPlusAB(NuComp,ReactComp)]

    # NOTE the variable `frozenReaction` is used to track which of the
    # potential reactions are frozen (i.e. only jumping not flowing).
    # This variable used to be called `discCompartment`. In a few
    # loops where there used to be a `compartIdx` iterator, these have
    # been changed to `reactionIdx` since these are looping over the
    # possible reactions not the compartments.
    frozenReaction = [False]*nRates
    for idx in range(nCompartments):
        for reactionIdx in range(nRates):
            if not EnforceDo[idx]:
                if DoDisc[idx] and compartInNu[reactionIdx][idx]:
                    frozenReaction[reactionIdx] = True
            else:
                if DoDisc[idx] and compartInNu[reactionIdx][idx]:
                    frozenReaction[reactionIdx] = True

    # initialise discrete sum compartments
    integralOfFiringTimes = [0.0]*nRates
    randTimes = [random.random() for _ in range(nRates)]


    tauArray = [Time(0.0)]*nRates

    overFlowAllocation = round(1000 * (t_max+dt)/dt + 1)

    # initialise solution arrays
    X = [[x0[i]] for i in range(nCompartments)]
    TauArr = [Time(0.0)]
    iters = 0

    # Track Absolute time
    AbsT = Time(0.0)
    ContT = Time(0.0)

    Xprev = x0
    Xcurr = x0

    # The `newlyDiscCompIndex` variable is used to track which
    # compartment has recently switched to being discrete. It was
    # previously called NewDiscCompartmemt, not to be confused with
    # NewfrozenReaction or frozenReaction which is a list of
    # booleans describing which reactions are frozen.
    newlyDiscCompIndex = None
    correctInteger = 0

    # import pdb; pdb.set_trace()
    while ContT < t_max:

        Dtau = dt
        Xprev = SystemState([x[iters] for x in X])
        Props = rates(Xprev, ContT)

        # Perform the Forward Euler Step
        dXdt = ComputedXdt(Props, nu, frozenReaction, nCompartments)

        # check if any states change in this step
        Dtau, correctInteger, DoDisc, frozenReaction, NewDoDisc, NewfrozenReaction, newlyDiscCompIndex = UpdateCompartmentRegime(dt, Xprev, Dtau, dXdt, Props, nu, SwitchingThreshold, DoDisc, EnforceDo, frozenReaction, compartInNu, nCompartments,nRates)

        # Only apply the forward Euler step if the compartment is continuous.
        Xcurr = SystemState([CompartmentValue(X[i][iters] + (0 if DoDisc[i] else Dtau * dXdt[i])) for i in range(nCompartments)])

        # Update the discrete compartments, if a state has just become discrete
        OriginalDoDisc = DoDisc[:]
        if correctInteger == 1:
            NewDoDisc, NewfrozenReaction = IsDiscrete(Xprev, SwitchingThreshold, DoDisc, EnforceDo, frozenReaction, compartInNu, nCompartments,nRates)
            frozenReaction = NewfrozenReaction[:]
            DoDisc = NewDoDisc[:]

        # Perform the Stochastic Loop
        stayWhile = any(DoDisc)

        AbsT = ContT
        DtauContStep = Dtau
        TimePassed = Time(0.0)

        firstStayWhileLoop = True
        while stayWhile:

            firstStayWhileLoop = False
            if TimePassed > 0:
                Props = rates(Xcurr, AbsT)

            integralStep = ComputeIntegralOfFiringTimes(Dtau, Props, rates, Xprev, Xcurr, AbsT)
            integralOfFiringTimes = ArrayPlusAB(integralOfFiringTimes,ArrayMultiplyAB(integralStep, frozenReaction))

            # If any of the components have just become discrete, we need to update the integralOfFiringTimes and randTimes
            if correctInteger == 1:
                for ii in range(nCompartments):
                    if newlyDiscCompIndex==ii and not EnforceDo[ii]:
                        for jj in range(nRates):
                            if compartInNu[jj][ii]:
                                frozenReaction[jj] = True
                                integralOfFiringTimes[jj] = 0.0
                                randTimes[jj] = random.random()

            firedReactions = [
                        ((0 > (rand - (1 - math.exp(-integral)))) and disc)
                        for rand, integral, disc in zip(randTimes, integralOfFiringTimes, frozenReaction)
                    ]

            if any(firedReactions):
                # Identify which reactions have fired
                tauArray = ComputeFiringTimes(firedReactions,integralOfFiringTimes,randTimes,Props,Dtau,nRates,integralStep)

                if num_non_zero(tauArray) > 0:

                    # Update the discrete compartments
                    Xcurr, Xprev, integralOfFiringTimes, integralStep, randTimes, TimePassed, AbsT, DtauMin = ImplementFiredReaction(tauArray ,integralOfFiringTimes,randTimes,Props,rates,integralStep,TimePassed, AbsT, X, iters, nu, dXdt, OriginalDoDisc, frozenReaction)

                    iters = iters + 1
                    for i in range(nCompartments):
                        X[i].append(Xcurr[i])

                    TauArr.append(AbsT)

                    Dtau = Dtau - DtauMin

                else:
                    stayWhile = False
            else:
                stayWhile = False

            if TimePassed >= DtauContStep:
                stayWhile = False

        iters = iters + 1
        ContT = Time(ContT + DtauContStep)

        TauArr.append(ContT)
        for i in range(len(X)):
            X[i].append(CompartmentValue(X[i][iters - 1] + (0 if DoDisc[i] else (DtauContStep - TimePassed) * dXdt[i])))

        if correctInteger == 1:
            pos = newlyDiscCompIndex
            X[pos][iters] = round(X[pos][iters])

            for jj in range(nRates):
                if compartInNu[jj][pos]:
                    frozenReaction[jj] = True
                    integralOfFiringTimes[jj] = 0.0
                    randTimes[jj] = random.random()
            newlyDiscCompIndex = None

    return Trajectory((X, TauArr))

def ComputeFiringTimes(firedReactions,integralOfFiringTimes,randTimes,Props,dt,nRates,integralStep):
    tauArray = [0.0] * nRates
    for kk in range(nRates):
        if firedReactions[kk]:
            Integral_t0_ti = integralStep[kk] - integralOfFiringTimes[kk]
            Integral = Integral_t0_ti - math.log((1 - randTimes[kk]))
            tauArray[kk] = Integral / Props[kk]

    return tauArray

def ImplementFiredReaction(tauArray ,integralOfFiringTimes,randTimes,Props,rates,integralStep,TimePassed, AbsT, X, iters, nu, dXdt, OriginalDoDisc,frozenReaction):
    tauArray = [float('inf') if tau == 0.0 else tau for tau in tauArray]

    DtauMin = min(tauArray)
    pos = tauArray.index(DtauMin)

    TimePassed = TimePassed + DtauMin
    AbsT = AbsT + DtauMin

    nCompartments = len(X)
    Xcurr = [X[i][iters] + nu[pos][i] + (0 if OriginalDoDisc[i] else DtauMin * dXdt[i]) for i in range(nCompartments)]
    Xprev = [X[i][iters] for i in range(nCompartments)]

    integralOfFiringTimes = [integral - step*disc for integral, step, disc in zip(integralOfFiringTimes, integralStep, frozenReaction)]
    integralStep = ComputeIntegralOfFiringTimes(DtauMin, Props, rates, Xprev, Xcurr, AbsT)
    integralOfFiringTimes = [integral + step*disc for integral, step, disc in zip(integralOfFiringTimes, integralStep, frozenReaction)]

    integralOfFiringTimes[pos] = 0.0
    randTimes[pos] = random.random()

    return Xcurr, Xprev, integralOfFiringTimes, integralStep, randTimes, TimePassed, AbsT, DtauMin

def ComputeIntegralOfFiringTimes(Dtau, Props, rates, Xprev, Xcurr, AbsT):
    # Integrate the cumulative wait times using trapezoid method
    integralStep = [Dtau * 0.5 * (p + r) for p, r in zip(Props, rates(Xcurr, AbsT + Dtau))]
    return integralStep


def ComputedXdt(
        Props: List[float],
        nu: List[List[float]],
        frozenReaction: List[bool],
        nCompartments: int) -> List[float]:
    """
    Compute the derivative of the state vector at time t for the
    continuous compartments by summing the contributions from each
    reaction.

    Args:
        Props: A list of the propensities of each reaction.
        nu: The stoichiometry matrix.
        frozenReaction: A list of booleans indicating which reactions
            are frozen.
        nCompartments: The number of compartments in the system.
    """
    return [sum(0 if frozenReaction[i] else Props[i] * nu[i][j]
                for i in range(len(Props)))
            for j in range(nCompartments)]


def UpdateCompartmentRegime(dt, Xprev, Dtau, dXdt, Props, nu, SwitchingThreshold, DoDisc, EnforceDo, frozenReaction, compartInNu, nCompartments,nRates):
    # check if any states change in this step
    NewDoDisc, NewfrozenReaction = IsDiscrete(Xprev, SwitchingThreshold, DoDisc, EnforceDo, frozenReaction, compartInNu, nCompartments,nRates)

    correctInteger = 0
    newlyDiscCompIndex = None
    x_step = [( 0 if isDisc else Dtau*dxi ) for dxi, isDisc in zip(dXdt, DoDisc)]
    if any([( (x+dxi  <= thresh) and  (not isDisc)) for x, isDisc, thresh, dxi in zip(Xprev,DoDisc,SwitchingThreshold,x_step)]):
        # Identify which compartment has just switched
        pos = 0
        possible_Dtau = [dt]
        for i, (x, isDisc, thresh, dxi) in enumerate(zip(Xprev,DoDisc,SwitchingThreshold, x_step)):
            if ((not isDisc) and  (x+dxi  <= thresh)):

                Xprev_pos = Xprev[i]
                rounded_Xprev_pos = math.ceil(Xprev_pos+x_step[i])
                dXdt_pos = dXdt[i]

                possible_Dtau.append(abs((rounded_Xprev_pos - Xprev_pos) / dXdt_pos))

        if len(possible_Dtau) > 1:
            Dtau = min(possible_Dtau)
            pos = possible_Dtau.index(Dtau)

            if pos > 0:
                    pos = pos - 1
                    # update the discrete compartment that has just switched
                    newlyDiscCompIndex = pos
                    correctInteger = 1
        else:
            frozenReaction = NewfrozenReaction
            DoDisc = NewDoDisc
            # correctInteer = 1
    else:
        frozenReaction = NewfrozenReaction
        DoDisc = NewDoDisc

    return Dtau, correctInteger, DoDisc, frozenReaction, NewDoDisc, NewfrozenReaction, newlyDiscCompIndex

def IsDiscrete(X, SwitchingThreshold, DoDisc, EnforceDo, frozenReaction, compartInNu, nCompartments, nRates):
    # Check if any compartments should be switched to continuous
    DoDiscTmp = [x <= threshold for x, threshold in zip(X, SwitchingThreshold)]

    for idx, x in enumerate(EnforceDo):
        if x == 1:
            DoDiscTmp[idx] = DoDisc[idx]

    are_equal = [DoDiscTmp[idx] == DoDisc[idx] for idx in range(nCompartments)]

    if all(are_equal):
        frozenReactionTmp = frozenReaction.copy()
    else:
        frozenReactionTmp = [False] * nRates

        for idx in range(nCompartments):
            for reactionIdx in range(nRates):
                if not EnforceDo[idx]:
                    if DoDiscTmp[idx] and compartInNu[reactionIdx][idx] == 1:
                        frozenReactionTmp[reactionIdx] = True
                else:
                    if DoDisc[idx] and compartInNu[reactionIdx][idx] == 1:
                        frozenReactionTmp[reactionIdx] = True


    return DoDiscTmp, frozenReactionTmp

# these are helper functions to make it readable
def ArraySubtractAB(ArrayA, ArrayB):
    AMinusB = [a - b for a, b in zip(ArrayA, ArrayB)]
    return AMinusB
def ArrayPlusAB(ArrayA, ArrayB):
    APlusB = [a + b for a, b in zip(ArrayA, ArrayB)]
    return APlusB
def MatrixSubtractAB(MatrixA,MatrixB):
    AMinusB = [[a - b for a, b in zip(row1, row2)] for row1, row2 in zip(MatrixA, MatrixB)]
    return AMinusB
def MatrixPlusAB(MatrixA,MatrixB):
    APlusB = [[a + b for a, b in zip(row1, row2)] for row1, row2 in zip(MatrixA, MatrixB)]
    return APlusB

def num_non_zero(array: Union[List[float], List[int], List[Time]]) -> int:
    return sum(1 for element in array if element != 0)

def MatrixDOTArray(Matrix,Array):
    result = [sum(row[i] * Array[i] for i in range(len(Array))) for row in Matrix]
    return result
def ArrayMultiplyAB(ArrayA, ArrayB):
    ArrayAB = [a * b for a, b in zip(ArrayA, ArrayB)]
    return ArrayAB
