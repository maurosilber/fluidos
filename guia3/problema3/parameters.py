import os
from pathlib import Path


def read_parameters(file: str | Path) -> str:
    with open(file) as f:
        return f.read()


def write_parameters(param_dir: Path, parameters: str):
    param_dir.mkdir(exist_ok=True)
    file = param_dir / "parameter.inp"

    if not read_parameters(file) == parameters:
        with open(file, "w") as f:
            f.write(parameters)

    os.link(solver_path, param_dir / "SOLVER")


if __name__ == "__main__":
    output_dir = Path("output")
    solver_path = Path("BOUSS")
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
