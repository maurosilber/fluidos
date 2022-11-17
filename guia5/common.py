from pathlib import Path

import numpy as np
import pandas as pd
import xarray

OUTPUT = dict(
    balance=["time", "<v^2>+<b^2>", "<omega^2>", "<j^2>"],
    energy=["time", "<v^2>", "<b^2>"],
    helicity=["time", "kinetic helicity", "magnetic helicity"],
    cross=["time", "<v.b>", "<a^2>"],
    kspectrum=["k", "Ev"],
    mspectrum=["k", "Eb"],
)


def load_global_output(path, name: str | list[str]):
    path = Path(path)
    if not isinstance(name, str):
        return pd.concat(
            (load_global_output(path, n) for n in name),
            axis=1,
            join="inner",
        )
    data = np.loadtxt(path / f"{name}.txt")
    columns = OUTPUT[name]
    return pd.DataFrame(data, columns=columns).set_index(columns[0])


def load_output(path, name: str | list[str], *, dt=1):
    path = Path(path)
    if not isinstance(name, str):
        return xarray.merge(load_output(path, n, dt=dt) for n in name)

    files = sorted(path.glob(f"{name}.*.txt"))
    index_name, *columns = OUTPUT[name]
    data = xarray.DataArray(
        [np.loadtxt(file)[:, 1:] for file in files],
        dims=["time", index_name, "data"],
        coords={
            "data": columns,
            "time": dt * np.arange(len(files)),
            index_name: np.loadtxt(files[0], usecols=0),
        },
    ).to_dataset("data")
    if len(data) == 1:
        name = next(iter(data.data_vars))
        data = data.to_array(name=name).squeeze(drop=True)
    return data
