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
        den_fsavg: str = 'ico128',
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
    den_fsavg: str, optional
        fsavg standard density. Choices are 'ico128', 'ico16', 'ico32', 'ico64'. Default 'ico128'.
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

    den_ind_fsavg = DENSITIES.index(den_fsavg)
    atlas_fsavg = list(ATLAS.items())[den_ind_fsavg]

    # set up transforms
    for hemi, h in zip(['lh', 'rh'], ['L', 'R']):

        #Set up some file paths
        fsavg_std_sphere_path = Path(HCPPIPEDIR, "global", "templates", "standard_mesh_atlases", "resample_fsaverage", f"{atlas_fsavg[0]}_std_sphere.{h}.{atlas_fsavg[1]}_fsavg_{h}.surf.gii")
        tpl_onavg_sphere_path = Path(cache_dir, f"tpl-onavg_hemi-{h}_den-{atlas[1]}_sphere-via_{atlas_fsavg[0]}.surf.gii")

        # initial individual register to fsaverage from fs_LR
        cmd = [
            "wb_command",
            "-surface-sphere-project-unproject",
            Path(native_dir, f"{subject}.{h}.sphere.MSMSulc.native.surf.gii").resolve(),
            Path(HCPPIPEDIR, "global", "templates", "standard_mesh_atlases", f"fsaverage.{h}_LR.spherical_std.164k_fs_LR.surf.gii").resolve(),
            Path(HCPPIPEDIR, "global", "templates", "standard_mesh_atlases", "resample_fsaverage", f"fs_LR-deformed_to-fsaverage.{h}.sphere.164k_fs_LR.surf.gii").resolve(),
            Path(outdir, f"{subject}.{h}.sphere.fsaverage_164k.native.surf.gii").resolve()
        ]
        subprocess.check_output(cmd)

        # resample midthickness surfaces to fsaverage
        cmd = [
            "wb_command",
            "-surface-resample",
            Path(native_dir, f"{subject}.{h}.{surface}.native.surf.gii").resolve(),
            Path(outdir, f"{subject}.{h}.sphere.fsaverage_164k.native.surf.gii").resolve(),
           fsavg_std_sphere_path.resolve(),
            "BARYCENTRIC",
            Path(outdir, f"{subject}.{h}.{surface}.{atlas_fsavg[1]}_fsavg_{h}.surf.gii").resolve()
        ]
        subprocess.check_output(cmd)

        # now register onavg sphere to fsaverage
        onavg_sph = Sphere(*read_geometry(str(Path(cache_dir, f"onavg-{den}", "surf", f"{hemi}.sphere.reg").resolve())))
        fs_sph = Sphere.from_gifti(str(fsavg_std_sphere_path.resolve()))
        #mid = Surface.from_gifti('surfaces/100206.L.midthickness.native.surf.gii')

        ## register
        if not Path.exists(tpl_onavg_sphere_path.resolve()):
            fs_to_onavg_mat = fs_sph.barycentric(onavg_sph.coords)
            fs_onavg_coords = fs_sph.coords.T @ fs_to_onavg_mat
            fs_onavg_coords = fs_onavg_coords.T

            fs_onavg_sph = Sphere(fs_onavg_coords, onavg_sph.faces)
            fs_onavg_sph.to_gifti(str(tpl_onavg_sphere_path.resolve()))

        # perform final register
        cmd = [
            "wb_command",
            "-surface-resample",
            Path(outdir, f"{subject}.{h}.{surface}.{atlas_fsavg[1]}_fsavg_{h}.surf.gii").resolve(),
            fsavg_std_sphere_path.resolve(),
            tpl_onavg_sphere_path.resolve(),
            "BARYCENTRIC",
            Path(outdir, f"{subject}.{h}.{surface}.onavg-{den}_via_{atlas_fsavg[0]}.surf.gii").resolve()
        ]
        subprocess.check_output(cmd)

def native_to_onavg(
        hcp_dir: Path | str, 
        subject: str, 
        surface: str = 'midthickness',
        den: str = 'ico128', 
        cache_dir: Path | str = None):
    
    # check for downloaded files
    if cache_dir is None:
        cache_dir = Path(os.environ['HOME'], 'onavg-template', parents=True, exist_ok=True).resolve()
    
    #native_dir = Path(hcp_dir, f"{subject}", "T1w", 'Native').resolve()
    outdir = Path(hcp_dir, f"{subject}", "onavg").resolve()
    if not outdir.exists():
        outdir.mkdir()

    if not Path.exists(Path(cache_dir, f'onavg-{den}')):
        get_onavg(cache_dir)
    
    subj_cache_dir = Path(cache_dir, 'individual', f'{subject}').resolve()
    if not subj_cache_dir.exists():
        subj_cache_dir.mkdir(parents=True)

    recon_dir = Path(hcp_dir, f"{subject}", "T1w", f"{subject}", "surf")
    for hemi, h in zip(['lh', 'rh'], ['L', 'R']):
        # now register onavg sphere to fsaverage
        onavg_sph = Sphere(*read_geometry(str(Path(cache_dir, f"onavg-{den}", "surf", f"{hemi}.sphere.reg").resolve())))
        native_sph = Sphere(*read_geometry(str(Path(recon_dir, f"{hemi}.sphere").resolve())))
        #mid = Surface.from_gifti('surfaces/100206.L.midthickness.native.surf.gii')

        ## register
        native_to_onavg_mat = native_sph.barycentric(onavg_sph.coords)
        native_onavg_coords = native_sph.coords.T @ native_to_onavg_mat
        native_onavg_coords = native_onavg_coords.T

        native_onavg_sph = Sphere(native_onavg_coords, onavg_sph.faces)
        native_onavg_sph.to_gifti(str(Path(subj_cache_dir, f"{subject}.{h}.sphere.onavg-{den}.surf.gii").resolve()))

        # convert recon-all native sphere to gifti for wb_command
        cmd = [
            "mris_convert",
            Path(recon_dir, f"{hemi}.sphere").resolve(),
            Path(subj_cache_dir, f"{subject}.{h}.sphere.native.surf.gii").resolve()
        ]
        subprocess.check_output(cmd)

        if surface == 'midthickness':
            cmd = [
                "mris_convert",
                Path(recon_dir, f"{hemi}.pial").resolve(),
                Path(subj_cache_dir, f"{subject}.{h}.pial.native.surf.gii").resolve()
            ]
            subprocess.check_output(cmd)

            cmd = [
                "mris_convert",
                Path(recon_dir, f"{hemi}.white").resolve(),
                Path(subj_cache_dir, f"{subject}.{h}.white.native.surf.gii").resolve()
            ]
            subprocess.check_output(cmd)
            cmd = [
                "wb_command",
                "-surface-average",
                Path(subj_cache_dir, f"{subject}.{h}.{surface}.native.surf.gii").resolve(),
                "-surf",
                Path(subj_cache_dir, f"{subject}.{h}.white.native.surf.gii").resolve(),
                "-surf",
                Path(subj_cache_dir, f"{subject}.{h}.pial.native.surf.gii").resolve()
            ]
            subprocess.check_output(cmd)
        
        else:
            cmd = [
                "mris_convert",
                Path(recon_dir, f"{hemi}.{surface}").resolve(),
                Path(subj_cache_dir, f"{subject}.{h}.{surface}.native.surf.gii").resolve()
            ]
            subprocess.check_output(cmd)
            
        # convert desired onavg-ico density to gifti
        # cmd = [
        #     "mris_convert",
        #     Path(cache_dir, f"onavg-{den}", "surf", f"{hemi}.sphere")
        # ]

        # resample native midthickness to onavg using above sphere
        cmd = [
            "wb_command",
            "-surface-resample",
            Path(subj_cache_dir, f"{subject}.{h}.{surface}.native.surf.gii").resolve(),
            Path(subj_cache_dir, f"{subject}.{h}.sphere.native.surf.gii").resolve(),
            Path(subj_cache_dir, f"{subject}.{h}.sphere.onavg-{den}.surf.gii").resolve(),
            "BARYCENTRIC",
            Path(outdir, f"{subject}.{h}.{surface}.onavg-{den}.v2.surf.gii").resolve()
        ]
        subprocess.check_output(cmd)

def nativeMNI_to_onavg(
        hcp_dir: Path | str, 
        subject: str, 
        surface: str = 'midthickness',
        den: str = 'ico128', 
        cache_dir: Path | str = None):
    
    hcp_dir = Path(hcp_dir).resolve()
    # check for downloaded files
    if cache_dir is None:
        cache_dir = Path(os.environ['HOME'], 'onavg-template', parents=True, exist_ok=True).resolve()
    if not Path.exists(Path(cache_dir, f'onavg-{den}')):
        get_onavg(cache_dir)

    recon_dir = Path(hcp_dir, f"{subject}", "T1w", f"{subject}", "surf")
    subj_cache_dir = Path(cache_dir, "individual", f"{subject}").resolve()
    if not subj_cache_dir.exists():
        subj_cache_dir.mkdir(parents=True)
    outdir = Path(hcp_dir, f"{subject}", "MNINonLinear", f"onavg").resolve()
    if not outdir.exists():
        outdir.mkdir()
    
    # check atlas
    den_ind = DENSITIES.index(den)
    atlas = list(ATLAS.items())[den_ind]

    for hemi, h in zip(['lh', 'rh'], ['L', 'R']):

        # register native to fsaverage
        cmd = [
            "mris_convert",
            Path(recon_dir, f"{hemi}.sphere.reg"),
            Path(subj_cache_dir, f"{subject}.{h}.sphere.reg.native.surf.gii"),
        ]
        subprocess.check_output(cmd)

        tplmain_dir = Path(cache_dir.parent,'tpl-onavg-main').resolve()
        trg_sphere = Path(tplmain_dir,f"tpl-onavg_hemi-{h}_den-{atlas[1]}_sphere.surf.gii") 
        cmd = [
            "wb_command",
            "-surface-resample",
            Path(hcp_dir, f"{subject}", "MNINonLinear", "Native", f"{subject}.{h}.{surface}.native.surf.gii"),
            Path(subj_cache_dir, f"{subject}.{h}.sphere.reg.native.surf.gii"),
            trg_sphere,
            "BARYCENTRIC",
            Path(outdir, f"{subject}.{h}.{surface}.onavg-{den}_v5.surf.gii"),
        ]
        subprocess.check_output(cmd)

def nativeMNI_to_fsaverage(
        hcp_dir: Path | str, 
        subject: str, 
        surface: str = 'midthickness',
        den: str = 'fsaverage5', 
        cache_dir: Path | str = None):
    
    # check atlas
    den_fsavg = ATLAS[den]

    hcp_dir = Path(hcp_dir).resolve()
    # check for downloaded files
    if cache_dir is None:
        cache_dir = Path(os.environ['HOME'], 'onavg-template', parents=True, exist_ok=True).resolve()

    recon_dir = Path(hcp_dir, f"{subject}", "T1w", f"{subject}", "surf")
    subj_cache_dir = Path(cache_dir, "individual", f"{subject}").resolve()
    FS_HOME = Path(os.environ['FREESURFER_HOME']).resolve()
    temp_dir = Path(FS_HOME, "subjects", f"{den}", "surf").resolve()

    outdir = Path(hcp_dir, f"{subject}", "MNINonLinear", f"{den}_{den_fsavg}").resolve()
    if not outdir.exists():
        outdir.mkdir()

    for hemi, h in zip(['lh', 'rh'], ['L', 'R']):
        # register native to fsaverage
        cmd = [
            "mris_convert",
            Path(recon_dir, f"{hemi}.sphere.reg"),
            Path(subj_cache_dir, f"{subject}.{h}.sphere.reg.native.surf.gii"),
        ]
        subprocess.check_output(cmd)

        cmd = [
            "mris_convert",
            Path(temp_dir, f"{hemi}.sphere"),
            Path(cache_dir, f"{den}.{h}.sphere.surf.gii"),
        ]
        subprocess.check_output(cmd)

        cmd = [
            "wb_command",
            "-surface-resample",
            Path(hcp_dir, f"{subject}", "MNINonLinear", "Native", f"{subject}.{h}.{surface}.native.surf.gii"),
            Path(subj_cache_dir, f"{subject}.{h}.sphere.reg.native.surf.gii"),
            Path(cache_dir, f"{den}.{h}.sphere.surf.gii"),
            "BARYCENTRIC",
            Path(outdir, f"{subject}.{h}.{surface}.{den}_{den_fsavg}.surf.gii"),
        ]
        subprocess.check_output(cmd)


def main():
    parser = ArgumentParser()
    parser.add_argument("HCPdir", help="path to post-HCP subject folders")
    parser.add_argument("subject", help="subject to register individual surface to onavg")
    parser.add_argument("surface", default="midthickness", metavar="surface", help="surface to register [ midthickness | pial | white ]")
    parser.add_argument("density", default='ico128', metavar="density", help="specify which onavg density to register individual surface [ ico128 | ico16 | ico32 | ico64]")
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
    print(f"registration success, saved to {hcp_dir}/{subject}/MNINonLinear/onavg/{subject}.L|R.{surface}.onavg-{den}.surf.gii")
    exit(0)
    
if __name__ == '__main__':
    main()