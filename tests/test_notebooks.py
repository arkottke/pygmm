"""Execute the example notebooks."""

import pathlib

import nbformat
import pytest
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

root = pathlib.Path(__file__).parent / "../examples"
fpaths = sorted(root.glob("*.ipynb"))


def idfn(val):
    if isinstance(val, pathlib.Path):
        return val.name


@pytest.mark.parametrize("fpath", fpaths, ids=idfn)
def test_notebook(fpath):
    """Execute each notebook, without modifying it."""
    # "python3" is the kernel that ipykernel provides for the interpreter running
    # the tests, so it does not depend on the kernels installed by the user or on
    # the kernelspec stored in the notebook. The executed copy is discarded.
    # https://nbclient.readthedocs.io/en/latest/client.html
    with fpath.open(encoding="utf-8") as fp:
        nb = nbformat.read(fp, as_version=4)

    client = NotebookClient(
        nb,
        timeout=600,
        kernel_name="python3",
        resources={"metadata": {"path": str(fpath.parent)}},
    )
    try:
        client.execute()
    except CellExecutionError as err:
        pytest.fail(f'Error executing the notebook "{fpath.name}":\n{err}')
