
# Description

The `jsf` package does not yet have an exact simulator available. This should be provided, even if the implementation is slow. It should strive to be as clear and obviously correct as possible. It must be quick enough to run all test cases though. It must handle the annoying edge case gracefully.

# Suggested solution

The new functionality should be provided in a module `jsf.exact`. The implemented function should have the following signature:

```python
def JumpSwitchFlowExact(
        x0: SystemState,
        rates: Callable[[SystemState, Time], List[float]],
        stoich: Dict[str, Any],
        t_max: Time,
        options: Dict[str, Any]) -> Tuple[List[List[CompartmentValue]], List[Time]]:
    pass
```

The `JumpSwitchFlowExact` function should coordinate a loop which builds up the process trajectory. The updating of the state at each step should be carried out by a separate function `_update` with the following signature:

```python
  JumpClock = NewType('JumpClock', float)
  ExtendedState = NewType('ExtendedState',
                          Tuple[SystemState, List[JumpClock], Time])

  def _update(
          x0: ExtendedState,
          delta_time: Time,
          rates: Callable[[SystemState, Time], List[float]],
          stoich: Dict[str, Any],
          options: Dict[str, Any]) -> ExtendedState:
      pass
```

The `_update` function should carry out a single step of forward-Euler. If this leads to a valid state of the system then that can be returned. Otherwise, we resolve this recursively: back-up to the time of the first event, apply the corresponding change in state, and then call `_update` again to resolve the remaining time until the end of the time step.

The two type of events that may occur are *jumps* and *switches*. There should be a function `_jump` that carries out the appropriate change to state upon that jump.

```python
  def _jump(
          x0: ExtendedState,
          reaction: int,
          stoich: Dict[str, Any],
          options: Dict[str, Any]) -> ExtendedState:
      pass
```

The resolution of the annoying edge case should be delegated to another function to make it easier to consider alternatives in the future.
