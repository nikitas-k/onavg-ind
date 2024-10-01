from . import Sphere
from .utils import read_geometry, get_onavg
import subprocess
from argparse import ArgumentParser
import requests
from pathlib import Path
import os

DENSITIES = [
    'ico128',
    'ico16',
    'ico32',
    'ico64',
]

ATLAS = {'fsaverage' : '164k',
    'fsaverage4' : '3k',
    'fsaverage5' : '10k',
    'fsaverage6' : '41k',
}

def register_to_onavg(
        hcp_dir: Path | str, 
        subject: str, 
        surface: str = 'midthickness',
        den: str = 'ico128', 
        cache_dir: Path | str = None):
    
    """
    Register individual surface <MNINonLinear/Native> to onavg. This function
    requires `wb_command` and assumes HCPpipelines `<https://github.com/Washington-University/HCPpipelines>`_
    to be cloned and correctly sourced (i.e., <HCPPIPEDIR> should be in your OS environment)

    Parameters
    ==========
    hcp_dir : pathlib.Path or str of path to file
        HCP-Pipelines post-processed subject directory. Expects the following directory and file structure
            hcp_dir---<subject>
                     |--MNINonLinear
                        |--Native
                           |--sub-001.<hemisphere>.sphere.MSMAll.native.surf.gii
                           |--sub-001.<hemisphere>.<surface>.surf.gii

    subject : str
        HCP-post processed subject folder
    surface : str, optional
        Surface to register to onavg standard. Default 'midthickness'.
    den : str, optional
        onavg standard density. Choices are 'ico128', 'ico16', 'ico32', 'ico64'. Default 'ico128'. 
    cache_dir : pathlib.Path or str, optional
        Directory to store onavg and standard freesurfer files for function. If None, defaults to ~/onavg-template. Default None.

    """
    if 'HCPPIPEDIR' not in os.environ:
        raise EnvironmentError('incorrectly sourced HCPPIPEDIR. check installation')
    HCPPIPEDIR = Path(os.environ['HCPPIPEDIR']).resolve()

    hcp_dir = Path(hcp_dir).resolve()
    # check for downloaded files
    if cache_dir is None:
        cache_dir = Path(os.environ['HOME'], 'onavg-template', parents=True, exist_ok=True).resolve()
    
    # download
    if not Path.exists(Path(cache_dir, f'onavg-{den}')):
        get_onavg(cache_dir)
    
    # check unzipped correctly
    if not Path.exists(Path(cache_dir, f'onavg-{den}')):
        raise RuntimeError(f'could not download files properly, check write permissions of folder: {cache_dir}')
    
    if not Path.exists(hcp_dir):
        raise ValueError('could not find HCPpipelines output directory, check')
    
    # check subject folder structure
    native_dir = Path(hcp_dir, subject, 'MNINonLinear', 'native', parents=False, exist_ok=False).resolve()
    outdir = Path(hcp_dir, subject, 'MNINonLinear', 'onavg', parents=True, exist_ok=True).resolve()
    if Path.exists(outdir) is False:
        Path.mkdir(outdir)

    # check atlas
    den_ind = DENSITIES.index(den)
    atlas = list(ATLAS.items())[den_ind]

    # set up transforms
    for hemi, h in zip(['lh', 'rh'], ['L', 'R']):
        # initial individual register to fsaverage from fs_LR
        cmd = [
            "wb_command",
            "-surface-sphere-project-unproject",
            Path(native_dir, f"{subject}.{h}.sphere.MSMAll.native.surf.gii").resolve(),
            Path(HCPPIPEDIR, "global", "templates", "standard_mesh_atlases", f"fsaverage.{h}_LR.spherical_std.{atlas[1]}_fs_LR.surf.gii").resolve(),
            Path(HCPPIPEDIR, "global", "templates", "standard_mesh_atlases", "resample_fsaverage", f"fs_LR-deformed_to-fsaverage.{h}.sphere.{atlas[1]}_fs_LR.surf.gii").resolve(),
            Path(outdir, f"{subject}.{h}.sphere.{atlas[0]}_{atlas[1]}.native.surf.gii").resolve()
        ]
        subprocess.check_output(cmd)

        # resample midthickness surfaces to fsaverage
        cmd = [
            "wb_command",
            "-surface-resample",
            Path(native_dir, f"{subject}.{h}.{surface}.native.surf.gii").resolve(),
            Path(outdir, f"{subject}.{h}.sphere.{atlas[0]}_{atlas[1]}.native.surf.gii").resolve(),
            Path(HCPPIPEDIR, "global", "templates", "standard_mesh_atlases", "resample_fsaverage", f"{atlas[0]}_std_sphere.{h}.{atlas[1]}_fsavg_{h}.surf.gii").resolve(),
            "BARYCENTRIC",
            Path(outdir, f"{subject}.{h}.{surface}.{atlas[1]}_fsavg_{h}.surf.gii").resolve()
        ]
        subprocess.check_output(cmd)

        # now register onavg sphere to fsaverage
        onavg_sph = Sphere(*read_geometry(str(Path(cache_dir, f"onavg-{den}", "surf", f"{hemi}.sphere.reg").resolve())))
        fs_sph = Sphere.from_gifti(str(Path(HCPPIPEDIR, "global", "templates", "standard_mesh_atlases", "resample_fsaverage" f"{atlas[0]}_std_sphere.{h}.{atlas[1]}_fsavg_{h}.surf.gii").resolve()))
        #mid = Surface.from_gifti('surfaces/100206.L.midthickness.native.surf.gii')

        ## register
        if not Path.exists(Path(cache_dir, f"tpl-onavg_hemi-{h}_den-{atlas[1]}_sphere-{atlas[0]}.surf.gii").resolve()):
            onavg_to_fs_mat = onavg_sph.barycentric(fs_sph.coords)
            onavg_fs_coords = onavg_sph.coords.T @ onavg_to_fs_mat
            onavg_fs_coords = onavg_fs_coords.T

            onavg_fs_sph = Sphere(onavg_fs_coords, fs_sph.faces)
            onavg_fs_sph.to_gifti(str(Path(cache_dir, f"tpl-onavg_hemi-{h}_den-{atlas[1]}_sphere-{atlas[0]}.surf.gii").resolve()))

        # perform final register
        cmd = [
            "wb_command",
            "-surface-resample",
            Path(outdir, f"{subject}.{h}.{surface}.{atlas[1]}_fsavg_{h}.surf.gii").resolve(),
            Path(HCPPIPEDIR, "global", "templates", "standard_mesh_atlases", "resample_fsaverage", f"{atlas[0]}_std_sphere.{h}.{atlas[1]}_fsavg_{h}.surf.gii").resolve(),
            Path(cache_dir, f"tpl-onavg_hemi-{h}_den-{atlas[1]}_sphere-{atlas[0]}.surf.gii").resolve(),
            "BARYCENTRIC",
            Path(outdir, f"{subject}.{h}.{surface}.onavg-{den}.surf.gii").resolve()
        ]
        subprocess.check_output(cmd)

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("HCPdir", help="path to post-HCP subject folders")
    parser.add_argument("subject", help="subject to register individual surface to onavg")
    parser.add_argument("surface", default="midthickness", metavar="midthickness pial white", help="surface to register")
    parser.add_argument("density", default='ico128', metavar="ico128 ico16 ico32 ico64", help="specify which onavg density to register individual surface")
    parser.add_argument("--cache-dir", default=None, help="download cache directory for onavg and HCP standard surfaces. defaults to ~/onavg-template")

    args = parser.parse_args()
    hcp_dir = args.HCPdir
    subject = args.subject
    surface = args.surface
    den = args.density
    if den not in DENSITIES:
        raise ValueError(f'error: unknown density {den}')
    
    cache_dir = args.cache_dir

    # main
    register_to_onavg(hcp_dir, subject, surface, den, cache_dir)
    print(f"registration success, saved to {hcp_dir}/{subject}/MNINonLinear/onavg/{subject}.?.{surface}.onavg-{den}.surf.gii")
    exit(0)