from __future__ import annotations

import os
from pathlib import Path


def read_parameters(file: str | Path) -> str:
    with open(file) as f:
        return f.read()


def write_parameters(param_dir: Path, parameters: str):
    param_dir.mkdir(exist_ok=True)
    file = param_dir / "parameter.inp"

    # Write parameters to param_dir if they:
    # - do not exist
    # - are different from current parameters
    try:
        assert read_parameters(file) == parameters
    except (AssertionError, FileNotFoundError):
        with open(file, "w") as f:
            f.write(parameters)
    
    # Hard-link solver to param_dir
    try:
        os.link("SOLVER", param_dir / "SOLVER")
    except FileExistsError:
        # TODO: que pasa si recompilo SOLVER?
        pass


if __name__ == "__main__":
    output_dir = Path("output")
    template = read_parameters("parameter.inp")

    def problema3(name, *, u1, N, gamma):
        parameters = template.format(
            u1=float(u1),
            bvfreq=float(N),
            gamma=float(gamma),
        )
        write_parameters(output_dir / name, parameters)

    problema3("punto_c", u1=0.1, N=2.0, gamma=10.0)

    for N in range(0, 7, 1):
        problema3(f"N{N}", u1=0.1, N=N, gamma=10)

    for gamma in range(5, 21, 5):
        problema3(f"gamma{gamma}", u1=0.1, N=2, gamma=gamma)
