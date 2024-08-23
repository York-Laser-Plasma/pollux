import numpy as np
import pathlib


def read_pollux(directory: str | pathlib.Path = "."):
    """Read a set of Pollux numpy files in `directory`

    Returns dictionary of time: data
    """

    directory = pathlib.Path(directory)

    result = {}
    for pollux_file in directory.glob("pollux_*.npz"):
        time = float(pollux_file.stem.split("_")[-1])
        with open(pollux_file, "rb") as f:
            result[time] = dict(np.load(f))

    # result keys will be in lexicographical order: reorder
    # numerically by time
    return {time: data for time, data in sorted(result.items())}
    
