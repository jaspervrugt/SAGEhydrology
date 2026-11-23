"""Explicit numerical failure; configuration and library errors stay fatal."""


class IntegrationFailure(ValueError):
    """A valid simulation request could not be integrated."""


def native_failure(message):
    numerical = (
        "Required timestep",
        "Required step",
        "Newton",
        "Singular",
        "Adaptive integration attempt",
        "Fixed-step attempt",
        "Nonfinite output",
        "Nonfinite sensitivity",
        "Nonfinite error",
        "Step limit",
        "BDF2 required",
        "RK2 cannot",
        "DP5 cannot",
    )
    return (
        IntegrationFailure(message)
        if message.startswith(numerical)
        else ValueError(message)
    )
