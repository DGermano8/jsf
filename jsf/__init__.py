import random
import math

def jsf(x0, rates, stoich, t_max, **kwargs):
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


def JumpSwitchFlowSimulator(x0, rates, stoich, t_max, options):
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
    # predefine and initialise the system
    # TODO - add default options

    # NOTE In the nu-matrix each row is a reaction and each column
    # describes the number items of that species used in the reaction.
    # There is one rate for each reaction and since there is a column
    # for each species we can use the length of the first column to
    # find the number of species in the system.
    nu = stoich["nu"]
    nRates = len(nu)
    nCompartments = len(nu[0])
    nuReactant = stoich["nuReactant"]

    dt = options["dt"]

    # This is to enable us to use a Boolean list here while still
    # accepting 0/1 as input. In future versions it would be nice to
    # require a Boolean list but that would break backwards
    # compatibility so we will leave it for now.
    EnforceDo = [(not (ed == 0)) for ed in options["EnforceDo"]]

    SwitchingThreshold = options["SwitchingThreshold"]

    DoDisc = [(x <= threshold and x==round(x)) for x, threshold in zip(x0, SwitchingThreshold)]


    # identify which compartment is in which reaction:
    NuComp = [[value != 0 for value in row] for row in nu]
    ReactComp = [[value != 0 for value in row] for row in nuReactant]
    compartInNu = [[value != 0 for value in row] for row in MatrixPlusAB(NuComp,ReactComp)]

    discCompartment = [False]*nRates
    for idx in range(nCompartments):
        for compartIdx in range(nRates):
            if not EnforceDo[idx]:
                if DoDisc[idx] and compartInNu[compartIdx][idx]:
                    discCompartment[compartIdx] = True
            else:
                if DoDisc[idx] and compartInNu[compartIdx][idx]:
                    discCompartment[compartIdx] = True

    # initialise discrete sum compartments
    integralOfFiringTimes = [0]*nRates
    randTimes = [random.random() for _ in range(nRates)]


    tauArray = [0]*nRates

    overFlowAllocation = round(1000 * (t_max+dt)/dt + 1)

    # initialise solution arrays
    X = [[x0[i]] for i in range(nCompartments)]
    TauArr = [0.0]
    iters = 0

    # Track Absolute time
    AbsT = 0
    ContT = 0

    Xprev = x0
    Xcurr = x0

    # NewDiscCompartmemt = [0.0] * nCompartments
    NewDiscCompartmemt = None
    correctInteger = 0

    # import pdb; pdb.set_trace()
    while ContT < t_max:

        Dtau = dt
        Xprev = [X[i][iters] for i in range(len(X))]
        Props = rates(Xprev, ContT)

        # Perform the Forward Euler Step
        dXdt = ComputedXdt(Xprev, Props, nu, discCompartment, nCompartments)

        # check if any states change in this step
        Dtau, correctInteger, DoDisc, discCompartment, NewDoDisc, NewdiscCompartment, NewDiscCompartmemt = UpdateCompartmentRegime(dt, Xprev, Dtau, dXdt, Props, nu, SwitchingThreshold, DoDisc, EnforceDo, discCompartment, compartInNu, nCompartments,nRates)

        # Only apply the forward Euler step if the compartment is continuous.
        Xcurr = [X[i][iters] + (0 if DoDisc[i] else Dtau * dXdt[i]) for i in range(nCompartments)]

        # Update the discrete compartments, if a state has just become discrete
        OriginalDoDisc = DoDisc[:]
        if correctInteger == 1:
            NewDoDisc, NewdiscCompartment = IsDiscrete(Xprev, SwitchingThreshold, DoDisc, EnforceDo, discCompartment, compartInNu, nCompartments,nRates)
            discCompartment = NewdiscCompartment[:]
            DoDisc = NewDoDisc[:]

        # Perform the Stochastic Loop
        stayWhile = any(DoDisc)

        AbsT = ContT
        DtauContStep = Dtau
        TimePassed = 0

        firstStayWhileLoop = True
        while stayWhile:

            firstStayWhileLoop = False
            if TimePassed > 0:
                Props = rates(Xcurr, AbsT)

            integralStep = ComputeIntegralOfFiringTimes(Dtau, Props, rates, Xprev, Xcurr, AbsT)
            integralOfFiringTimes = ArrayPlusAB(integralOfFiringTimes,ArrayMultiplyAB(integralStep, discCompartment))

            # If any of the components have just become discrete, we need to update the integralOfFiringTimes and randTimes
            if correctInteger == 1:
                for ii in range(nCompartments):
                    if NewDiscCompartmemt==ii and not EnforceDo[ii]:
                        for jj in range(nRates):
                            if compartInNu[jj][ii]:
                                discCompartment[jj] = 1
                                integralOfFiringTimes[jj] = 0.0
                                randTimes[jj] = random.random()

            firedReactions = [
                        (0 > rand - (1 - math.exp(-integral))) * disc
                        for rand, integral, disc in zip(randTimes, integralOfFiringTimes, discCompartment)
                    ]

            if NNZ(firedReactions) > 0:
                # Identify which reactions have fired
                tauArray = ComputeFiringTimes(firedReactions,integralOfFiringTimes,randTimes,Props,Dtau,nRates,integralStep)

                if NNZ(tauArray) > 0:

                    # Update the discrete compartments
                    Xcurr, Xprev, integralOfFiringTimes, integralStep, randTimes, TimePassed, AbsT, DtauMin = ImplementFiredReaction(tauArray ,integralOfFiringTimes,randTimes,Props,rates,integralStep,TimePassed, AbsT, X, iters, nu, dXdt, OriginalDoDisc, discCompartment)

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
        ContT = ContT + DtauContStep

        TauArr.append(ContT)
        for i in range(len(X)):
            X[i].append(X[i][iters - 1] + (0 if DoDisc[i] else (DtauContStep - TimePassed) * dXdt[i]))

        if correctInteger == 1:
            pos = NewDiscCompartmemt
            X[pos][iters] = round(X[pos][iters])

            for jj in range(nRates):
                if compartInNu[jj][pos]:
                    discCompartment[jj] = 1
                    integralOfFiringTimes[jj] = 0.0
                    randTimes[jj] = random.random()
            NewDiscCompartmemt = None

    return X, TauArr

def ComputeFiringTimes(firedReactions,integralOfFiringTimes,randTimes,Props,dt,nRates,integralStep):
    tauArray = [0.0] * nRates
    for kk in range(nRates):
        if firedReactions[kk]:
            Integral_t0_ti = integralStep[kk] - integralOfFiringTimes[kk]
            Integral = Integral_t0_ti - math.log((1 - randTimes[kk]))
            tauArray[kk] = Integral / Props[kk]

    return tauArray

def ImplementFiredReaction(tauArray ,integralOfFiringTimes,randTimes,Props,rates,integralStep,TimePassed, AbsT, X, iters, nu, dXdt, OriginalDoDisc,discCompartment):
    tauArray = [float('inf') if tau == 0.0 else tau for tau in tauArray]

    DtauMin = min(tauArray)
    pos = tauArray.index(DtauMin)

    TimePassed = TimePassed + DtauMin
    AbsT = AbsT + DtauMin

    nCompartments = len(X)
    Xcurr = [X[i][iters] + nu[pos][i] + (0 if OriginalDoDisc[i] else DtauMin * dXdt[i]) for i in range(nCompartments)]
    Xprev = [X[i][iters] for i in range(nCompartments)]

    integralOfFiringTimes = [integral - step*disc for integral, step, disc in zip(integralOfFiringTimes, integralStep, discCompartment)]
    integralStep = ComputeIntegralOfFiringTimes(DtauMin, Props, rates, Xprev, Xcurr, AbsT)
    integralOfFiringTimes = [integral + step*disc for integral, step, disc in zip(integralOfFiringTimes, integralStep, discCompartment)]

    integralOfFiringTimes[pos] = 0.0
    randTimes[pos] = random.random()

    return Xcurr, Xprev, integralOfFiringTimes, integralStep, randTimes, TimePassed, AbsT, DtauMin

def ComputeIntegralOfFiringTimes(Dtau, Props, rates, Xprev, Xcurr, AbsT):
    # Integrate the cumulative wait times using trapezoid method
    integralStep = [Dtau * 0.5 * (Props[i] + rates(Xcurr, AbsT + Dtau)[i]) for i in range(len(Props))]
    return integralStep


def ComputedXdt(Xprev, Props, nu, discCompartment, nCompartments):
    """
    Compute the derivative of the state vector at time t for the
    continuous compartments by summing the contributions from each
    reaction.
    """
    return [sum(0 if discCompartment[i] else Props[i] * nu[i][j] for i in range(len(Props))) for j in range(nCompartments)]


def UpdateCompartmentRegime(dt, Xprev, Dtau, dXdt, Props, nu, SwitchingThreshold, DoDisc, EnforceDo, discCompartment, compartInNu, nCompartments,nRates):
    # check if any states change in this step
    NewDoDisc, NewdiscCompartment = IsDiscrete(Xprev, SwitchingThreshold, DoDisc, EnforceDo, discCompartment, compartInNu, nCompartments,nRates)

    correctInteger = 0
    NewDiscCompartmemt = None
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
                    NewDiscCompartmemt = pos
                    correctInteger = 1
        else:
            discCompartment = NewdiscCompartment
            DoDisc = NewDoDisc
            # correctInteer = 1
    else:
        discCompartment = NewdiscCompartment
        DoDisc = NewDoDisc

    return Dtau, correctInteger, DoDisc, discCompartment, NewDoDisc, NewdiscCompartment, NewDiscCompartmemt

def IsDiscrete(X, SwitchingThreshold, DoDisc, EnforceDo, discCompartment, compartInNu, nCompartments, nRates):
    # Check if any compartments should be switched to continuous
    DoDiscTmp = [x <= threshold for x, threshold in zip(X, SwitchingThreshold)]

    for idx, x in enumerate(EnforceDo):
        if x == 1:
            DoDiscTmp[idx] = DoDisc[idx]

    are_equal = [DoDiscTmp[idx] == DoDisc[idx] for idx in range(nCompartments)]

    if all(are_equal):
        discCompartmentTmp = discCompartment.copy()
    else:
        discCompartmentTmp = [0] * nRates

        for idx in range(nCompartments):
            for compartIdx in range(nRates):
                if not EnforceDo[idx]:
                    if DoDiscTmp[idx] and compartInNu[compartIdx][idx] == 1:
                        discCompartmentTmp[compartIdx] = 1
                else:
                    if DoDisc[idx] and compartInNu[compartIdx][idx] == 1:
                        discCompartmentTmp[compartIdx] = 1


    return DoDiscTmp, discCompartmentTmp

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
def NNZ(Array):
    # NumberOfNonZeros
    non_zero_count = sum(1 for element in Array if element != 0)
    return non_zero_count
def MatrixDOTArray(Matrix,Array):
    result = [sum(row[i] * Array[i] for i in range(len(Array))) for row in Matrix]
    return result
def ArrayMultiplyAB(ArrayA, ArrayB):
    ArrayAB = [a * b for a, b in zip(ArrayA, ArrayB)]
    return ArrayAB
