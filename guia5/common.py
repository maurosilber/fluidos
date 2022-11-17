from pathlib import Path

import numpy as np
import pandas as pd
import xarray


class Output:
    def __init__(self, path, dt=1, **column_names: list[str]):
        self.path = Path(path)
        self.dt = dt
        self.column_names = column_names

    def load_global_output(self, *names: str):
        if len(names) > 1:
            return pd.concat(
                (self.load_global_output(n) for n in names),
                axis=1,
                join="inner",
            )

        name = names[0]
        data = np.loadtxt(self.path / f"{name}.txt")
        columns = self.column_names[name]
        return pd.DataFrame(data, columns=columns).set_index(columns[0])

    def load_output(self, *names: str):
        if len(names) > 1:
            return xarray.merge(self.load_output(n) for n in names)

        name = names[0]
        files = sorted(self.path.glob(f"{name}.*.txt"))
        index_name, *columns = self.column_names[name]
        data = xarray.DataArray(
            [np.loadtxt(file)[:, 1:] for file in files],
            dims=["time", index_name, "data"],
            coords={
                "data": columns,
                "time": self.dt * np.arange(len(files)),
                index_name: np.loadtxt(files[0], usecols=0),
            },
        ).to_dataset("data")
        if len(data) == 1:
            name = next(iter(data.data_vars))
            data = data.to_array(name=name).squeeze(drop=True)
        return data
