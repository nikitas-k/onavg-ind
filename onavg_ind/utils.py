from collections import OrderedDict
import warnings
import requests
import zipfile
import io
import os
import hashlib

import numpy as np

HASH = "87a17ac12f0634aca1f3ab615edc5a99"

def _fread3(fobj):
    """Read a 3-byte int from an open binary file object
    Parameters
    ----------
    fobj : file
        File descriptor
    Returns
    -------
    n : int
        A 3 byte int
    """
    b1, b2, b3 = np.fromfile(fobj, ">u1", 3)
    return (b1 << 16) + (b2 << 8) + b3

def _read_volume_info(fobj):
    """Helper for reading the footer from a surface file."""
    volume_info = OrderedDict()
    head = np.fromfile(fobj, ">i4", 1)
    if not np.array_equal(head, [20]):  # Read two bytes more
        head = np.concatenate([head, np.fromfile(fobj, ">i4", 2)])
        if not np.array_equal(head, [2, 0, 20]) and not np.array_equal(
            head, [2, 1, 20]
        ):
            warnings.warn("Unknown extension code.")
            return volume_info
        head = [2, 0, 20]

    volume_info["head"] = head
    for key in [
        "valid",
        "filename",
        "volume",
        "voxelsize",
        "xras",
        "yras",
        "zras",
        "cras",
    ]:
        pair = fobj.readline().decode("utf-8").split("=")
        if pair[0].strip() != key or len(pair) != 2:
            raise IOError("Error parsing volume info.")
        if key in ("valid", "filename"):
            volume_info[key] = pair[1].strip()
        elif key == "volume":
            volume_info[key] = np.array(pair[1].split()).astype(int)
        else:
            volume_info[key] = np.array(pair[1].split()).astype(float)
    # Ignore the rest
    return volume_info

def read_geometry(filepath, read_metadata=False, read_stamp=False):
    """Read a triangular format Freesurfer surface mesh.
    Parameters
    ----------
    filepath : str
        Path to surface file.
    read_metadata : bool, optional
        If True, read and return metadata as key-value pairs.
        Valid keys:
        * 'head' : array of int
        * 'valid' : str
        * 'filename' : str
        * 'volume' : array of int, shape (3,)
        * 'voxelsize' : array of float, shape (3,)
        * 'xras' : array of float, shape (3,)
        * 'yras' : array of float, shape (3,)
        * 'zras' : array of float, shape (3,)
        * 'cras' : array of float, shape (3,)
    read_stamp : bool, optional
        Return the comment from the file
    Returns
    -------
    coords : numpy array
        nvtx x 3 array of vertex (x, y, z) coordinates.
    faces : numpy array
        nfaces x 3 array of defining mesh triangles.
    volume_info : OrderedDict
        Returned only if `read_metadata` is True.  Key-value pairs found in the
        geometry file.
    create_stamp : str
        Returned only if `read_stamp` is True.  The comment added by the
        program that saved the file.
    """
    volume_info = OrderedDict()

    TRIANGLE_MAGIC = 16777214

    with open(filepath, "rb") as fobj:
        magic = _fread3(fobj)

        #if magic == TRIANGLE_MAGIC:  # Triangle file
        create_stamp = fobj.readline().rstrip(b"\n").decode("utf-8")
        test_dev = fobj.peek(1)[:1]
        if test_dev == b"\n":
            fobj.readline()
        vnum = np.fromfile(fobj, ">i4", 1)[0]
        fnum = np.fromfile(fobj, ">i4", 1)[0]
        coords = np.fromfile(fobj, ">f4", vnum * 3).reshape(vnum, 3)
        faces = np.fromfile(fobj, ">i4", fnum * 3).reshape(fnum, 3)

        if read_metadata:
            volume_info = _read_volume_info(fobj)
        # else:
        #     raise ValueError(
        #         f"{filepath} does not appear to be a Freesurfer surface (triangle file)"
        #     )

    coords = coords.astype(float)  # XXX: due to mayavi bug on mac 32bits

    ret = (coords, faces)
    if read_metadata:
        if len(volume_info) == 0:
            warnings.warn("No volume information contained in the file")
        ret += (volume_info,)
    if read_stamp:
        ret += (create_stamp,)

    return ret

def get_onavg(cache_dir):
    # get onavg templates
    url = "https://www.dropbox.com/s/rvg6ui3dqjx7v8i/tpl-onavg_freesurfer.zip?dl=1"
    response = requests.get(url)
    response.raise_for_status()
    
    with zipfile.ZipFile(io.BytesIO(response.content)) as z:
        z.extractall(cache_dir)