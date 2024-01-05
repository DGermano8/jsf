from typing import List, NewType

CompartmentValue = NewType('CompartmentValue', float)
SystemState = NewType('SystemState', List[CompartmentValue])
Time = NewType('Time', float)
