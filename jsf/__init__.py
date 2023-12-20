import random
import math
# Import the fibonnaci function as defined in ./cython/hello.pyx. The
# actual definition of the function is in ./src/hello.c and the header
# is ./src/hello.h.
from jsf.hello import py_fibonacci as fib


def JumpSwitchFlowSimulator(x0, rates, stoich, times, options):
    """
    x0 :: State                 # state at time zero.
    rates :: State -> Time -> [Rate] # rate function.
    stoich :: Structure
    times :: Number             # think "final time"
    options :: Structure
    """
    # predefine and initialise the system
    # TODO - add default options

    X0 = x0
    nu = stoich["nu"]
    nuReactant = stoich["nuReactant"]

    tFinal = times
    dt = options["dt"]
    EnforceDo = options["EnforceDo"]
    SwitchingThreshold = options["SwitchingThreshold"]

    DoDisc = stoich["DoDisc"]
    DoDisc = [(1 if x <= threshold and x==round(x) else 0) for x, threshold in zip(X0, SwitchingThreshold)]
    DoCont = ArraySubtractAB([1]*len(DoDisc), DoDisc)

    nRates = len(nu)
    nCompartments = len(nu[0])

    # identify which compartment is in which reaction:
    NuComp = [[value != 0 for value in row] for row in nu]
    ReactComp = [[value != 0 for value in row] for row in nuReactant]
    compartInNu = [[value != 0 for value in row] for row in MatrixPlusAB(NuComp,ReactComp)]

    discCompartment = [0]*nRates
    for idx in range(nCompartments):
        for compartIdx in range(nRates):
            if  EnforceDo[idx]==0:
                if DoDisc[idx]==1 and compartInNu[compartIdx][idx]:
                    discCompartment[compartIdx] = 1
            else:
                if DoDisc[idx]==1 and compartInNu[compartIdx][idx]:
                    discCompartment[compartIdx] = 1
    contCompartment = ArraySubtractAB([1]*nRates,discCompartment)
    # initialise discrete sum compartments
    integralOfFiringTimes = [0]*nRates
    randTimes = [random.random() for _ in range(nRates)]


    tauArray = [0]*nRates

    overFlowAllocation = round(1000 * (tFinal+dt)/dt + 1)

    # initialise solution arrays
    X = [[X0[i]] for i in range(nCompartments)]
    TauArr = [0.0]
    iters = 0

    # Track Absolute time
    AbsT = 0
    ContT = 0

    Xprev = X0
    Xcurr = X0

    # NewDiscCompartmemt = [0.0] * nCompartments
    NewDiscCompartmemt = None
    correctInteger = 0

    # import pdb; pdb.set_trace()
    while ContT < tFinal:

        Dtau = dt
        Xprev = [X[i][iters] for i in range(len(X))]
        Props = rates(Xprev, ContT)

        # Perform the Forward Euler Step
        dXdt = ComputedXdt(Xprev, Props, nu, contCompartment, nCompartments)

        # check if any states change in this step
        Dtau, correctInteger, DoDisc, DoCont, discCompartment, contCompartment, NewDoDisc, NewDoCont, NewdiscCompartment, NewcontCompartment, NewDiscCompartmemt = UpdateCompartmentRegime(dt, Xprev, Dtau, dXdt, Props, nu, SwitchingThreshold, DoDisc, DoCont, EnforceDo, discCompartment, contCompartment, compartInNu, nCompartments,nRates)

        Xcurr = [X[i][iters] + Dtau * (dXdt[i] * DoCont[i]) for i in range(len(X))]

        # Update the discrete compartments, if a state has just become discrete
        OriginalDoCont = DoCont[:]
        if correctInteger == 1:
            NewDoDisc, NewDoCont, NewdiscCompartment, NewcontCompartment = IsDiscrete(Xprev, SwitchingThreshold, DoDisc, DoCont, EnforceDo, discCompartment, contCompartment, compartInNu, nCompartments,nRates)
            contCompartment = NewcontCompartment[:]
            discCompartment = NewdiscCompartment[:]
            DoCont = NewDoCont[:]
            DoDisc = NewDoDisc[:]


        # Perform the Stochastic Loop
        stayWhile = True if NNZ(DoCont) != nCompartments else False

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
                    Xcurr, Xprev, integralOfFiringTimes, integralStep, randTimes, TimePassed, AbsT, DtauMin = ImplementFiredReaction(tauArray ,integralOfFiringTimes,randTimes,Props,rates,integralStep,TimePassed, AbsT, X, iters, nu, dXdt, OriginalDoCont, discCompartment)

                    iters = iters + 1
                    for i in range(len(X)):
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
            X[i].append(X[i][iters - 1] + (DtauContStep - TimePassed) * (dXdt[i] * DoCont[i]))

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

            Integral_t0_ti = (-1.0*(integralOfFiringTimes[kk] - integralStep[kk]))
            Integral = Integral_t0_ti - math.log((1 - randTimes[kk]))
            tau_val_1 = Integral / (Props[kk])

            tau_val = tau_val_1

            tauArray[kk] = tau_val

    return tauArray

def ImplementFiredReaction(tauArray ,integralOfFiringTimes,randTimes,Props,rates,integralStep,TimePassed, AbsT, X, iters, nu, dXdt, OriginalDoCont,discCompartment):
    tauArray = [float('inf') if tau == 0.0 else tau for tau in tauArray]

    DtauMin = min(tauArray)
    pos = tauArray.index(DtauMin)

    TimePassed = TimePassed + DtauMin
    AbsT = AbsT + DtauMin

    Xcurr = [X[i][iters] + nu[pos][i] + DtauMin * (dXdt[i] * OriginalDoCont[i]) for i in range(len(X))]
    Xprev = [X[i][iters] for i in range(len(X))]

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

def ComputedXdt(Xprev, Props, nu, contCompartment, nCompartments):
    dXdt = [sum(Props[i] * contCompartment[i] * nu[i][j] for i in range(len(Props))) for j in range(len(nu[0]))]

    return dXdt

def ORIGINAL_UpdateCompartmentRegime(dt, Xprev, Dtau, Props, nu, SwitchingThreshold, DoDisc, DoCont, EnforceDo, discCompartment, contCompartment, compartInNu, nCompartments,nRates):
    # check if any states change in this step
    NewDoDisc, NewDoCont, NewdiscCompartment, NewcontCompartment = IsDiscrete(Xprev, SwitchingThreshold, DoDisc, DoCont, EnforceDo, discCompartment, contCompartment, compartInNu, nCompartments,nRates)

    correctInteger = 0
    NewDiscCompartmemt = [0] * nCompartments

    # Identify if a state has just become discrete
    
    if any([((x > (thresh)) and  (x  <= (thresh+1)) and  isCont) for x, isCont, thresh in zip(Xprev,DoCont,SwitchingThreshold)]):
    # if any([((x > (thresh)) and  (x+Dtau*dxi*isCont  <= (thresh)) and  isCont) for x, isCont, thresh, dxi in zip(Xprev,DoCont,SwitchingThreshold,dXdt)]):
        # Identify which compartment has just switched
        pos = None
        for i, (x, isCont, thresh) in enumerate(zip(Xprev,DoCont,SwitchingThreshold)):
            if (x <= (thresh+1) and  isCont):
                pos = i
                break

        if pos is not None:
            # Here, we identify the time to the next integer
            dXdt = [sum(Props[i] * contCompartment[i] * nu[i][j] for i in range(len(Props))) for j in range(len(nu[0]))]

            Xprev_pos = Xprev[pos]
            rounded_Xprev_pos = math.floor(Xprev_pos)
            dXdt_pos = dXdt[pos]

            Dtau = min(dt, abs((rounded_Xprev_pos - Xprev_pos) / dXdt_pos))

            # If the time to the next integer is less than the time step, we need to move the mesh
            if Dtau < dt:
                NewDiscCompartmemt[pos] = 1
                correctInteger = 1
        else:
            contCompartment = NewcontCompartment
            discCompartment = NewdiscCompartment
            DoCont = NewDoCont
            DoDisc = NewDoDisc
            # correctInteer = 1
    else:
        contCompartment = NewcontCompartment
        discCompartment = NewdiscCompartment
        DoCont = NewDoCont
        DoDisc = NewDoDisc

    return Dtau, correctInteger, DoDisc, DoCont, discCompartment, contCompartment, NewDoDisc, NewDoCont, NewdiscCompartment, NewcontCompartment, NewDiscCompartmemt

def UpdateCompartmentRegime(dt, Xprev, Dtau, dXdt, Props, nu, SwitchingThreshold, DoDisc, DoCont, EnforceDo, discCompartment, contCompartment, compartInNu, nCompartments,nRates):
    # check if any states change in this step
    NewDoDisc, NewDoCont, NewdiscCompartment, NewcontCompartment = IsDiscrete(Xprev, SwitchingThreshold, DoDisc, DoCont, EnforceDo, discCompartment, contCompartment, compartInNu, nCompartments,nRates)

    correctInteger = 0
    NewDiscCompartmemt = None
    x_step = [Dtau*dxi*isCont for dxi, isCont in zip(dXdt, DoCont)]
    if any([( (x+dxi  <= (thresh)) and  isCont==1) for x, isCont, thresh, dxi in zip(Xprev,DoCont,SwitchingThreshold,x_step)]):
        # Identify which compartment has just switched
        pos = 0
        possible_Dtau = [dt]
        for i, (x, isCont, thresh, dxi) in enumerate(zip(Xprev,DoCont,SwitchingThreshold, x_step)):
            # if (x > (thresh) and  isCont==1 and  (x+dxi  <= (thresh)) ):
            if (isCont==1 and  (x+dxi  <= (thresh)) ):
                # print('pos = ', i)
                # print('x = ', x)

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
            contCompartment = NewcontCompartment
            discCompartment = NewdiscCompartment
            DoCont = NewDoCont
            DoDisc = NewDoDisc
            # correctInteer = 1
    else:
        contCompartment = NewcontCompartment
        discCompartment = NewdiscCompartment
        DoCont = NewDoCont
        DoDisc = NewDoDisc

    return Dtau, correctInteger, DoDisc, DoCont, discCompartment, contCompartment, NewDoDisc, NewDoCont, NewdiscCompartment, NewcontCompartment, NewDiscCompartmemt

def IsDiscrete(X, SwitchingThreshold, DoDisc, DoCont, EnforceDo, discCompartment, contCompartment, compartInNu, nCompartments, nRates):
    # Check if any compartments should be switched to continuous
    DoDiscTmp = [(1 if x <= threshold else 0) for x, threshold in zip(X, SwitchingThreshold)]
    DoContTmp = [1 - x for x in DoDiscTmp]

    for idx, x in enumerate(EnforceDo):
        if x == 1:
            DoDiscTmp[idx] = DoDisc[idx]
            DoContTmp[idx] = DoCont[idx]

    are_equal = [DoDiscTmp[idx] == DoDisc[idx] for idx in range(nCompartments)]

    if sum(are_equal) == nCompartments:
        discCompartmentTmp = discCompartment.copy()
        contCompartmentTmp = contCompartment.copy()
    else:
        discCompartmentTmp = [0] * nRates

        for idx in range(nCompartments):
            for compartIdx in range(nRates):
                if EnforceDo[idx] == 0:
                    if DoDiscTmp[idx] == 1 and compartInNu[compartIdx][idx] == 1:
                        discCompartmentTmp[compartIdx] = 1
                else:
                    if DoDisc[idx] == 1 and compartInNu[compartIdx][idx] == 1:
                        discCompartmentTmp[compartIdx] = 1

    contCompartmentTmp = [1 - x for x in discCompartmentTmp]


    return DoDiscTmp, DoContTmp, discCompartmentTmp, contCompartmentTmp

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
