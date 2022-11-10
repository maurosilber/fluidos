from pathlib import Path

output_dir = Path("output")

with open("parameter.inp") as f:
    parameters = f.read()

for N in range(1, 15, 3):
    out_file = output_dir / f"N_{N}/parameter.inp"
    out_file.parent.mkdir(exist_ok=True)
    with open(out_file, "w") as f:
        f.write(parameters.format(N=float(N)))
