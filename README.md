# Pollux

Pollux hydro-code

There is currently no licence: only members of the York Plasma
Institute may use or copy this code.

## Building

You'll need:

- Fortran compiler (like `gfortran`)
- CMake
- `zip` (for numpy `npz` output)

```bash
$ cmake . -B build
$ cmake --build build
```

This creates `build/pollux`, the executable

## Reading numpy data

Numpy output is one file per timestep. Use the `read_pollux` function
to read the data into a single `dict`:

```python
from read_pollux import read_pollux

data = read_pollux()
print(data.keys())
# dict_keys([1e-13, 4.746e-13, 1.812e-12, ...

print(data[1e-13].keys())
# dict_keys(['rho', 'Te', 'Ti', 'u', 'v'])
```
