"""
    The possible exit codes of the pipeline components.
"""
from enum import Enum, unique


@unique
class ExitCode(Enum):
    Success = 0

    MaxsatModelUnsat = 10
    NoAbstractionUnderComplexityBound = 11

    AbstractionFailsOnTestInstances = 20

    FeatureGenerationUnknownError = 30

    OutOfMemory = 100
    OutOfTime = 101
