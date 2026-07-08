"""
Entry-point shims for pip console_scripts.
Each function loads the corresponding extensionless script from pex_scripts/
via importlib so the scripts keep their clean names on the CLI.
"""
import importlib.util
import importlib.machinery
import os as _os
import sys as _sys


def _run(name):
    # Add TEMPy source dir so bare imports like 'from Vector import *' work.
    # Use the known path relative to this file rather than relying on import.
    _repo_root = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
    _tempy_src = _os.path.join(_repo_root, 'TEMPy_py3', 'src', 'TEMPy')
    if _os.path.isdir(_tempy_src) and _tempy_src not in _sys.path:
        _sys.path.insert(0, _tempy_src)

    path = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'pex_scripts', name)
    loader = importlib.machinery.SourceFileLoader(name, path)
    spec = importlib.util.spec_from_file_location(name, path, loader=loader)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.main()


def etomo_from_aretomo2():            _run('etomo_from_aretomo2')
def make_fid_out_of_3d_model():       _run('make_fid_out_of_3d_model')
def pex_add_maps():                   _run('pex_add_maps')
def pex_angular_distribution():       _run('pex_angular_distribution')
def pex_bin_model_files():            _run('pex_bin_model_files')
def pex_clean_by_planes():            _run('pex_clean_by_planes')
def pex_clean_by_tilt_ang():          _run('pex_clean_by_tilt_ang')
def pex_clean_using_chim_markers():   _run('pex_clean_using_chim_markers')
def pex_clean_using_imod_model():     _run('pex_clean_using_imod_model')
def pex_combine_model_files():        _run('pex_combine_model_files')
def pex_combine_prm_files():          _run('pex_combine_prm_files')
def pex_get_cccs():                   _run('pex_get_cccs')
def pex_get_hexons():                 _run('pex_get_hexons')
def pex_get_icos_symm_pcles():        _run('pex_get_icos_symm_pcles')
def pex_get_pentons():                _run('pex_get_pentons')
def pex_get_tetra():                  _run('pex_get_tetra')
def pex_get_zeros_of_ctf():           _run('pex_get_zeros_of_ctf')
def pex_linear_pick():                _run('pex_linear_pick')
def pex_make_bilds():                 _run('pex_make_bilds')
def pex_make_chim_markers():          _run('pex_make_chim_markers')
def pex_make_chim_markers_from_prm(): _run('pex_make_chim_markers_from_prm')
def pex_make_cylinder():              _run('pex_make_cylinder')
def pex_make_link_markers():          _run('pex_make_link_markers')
def pex_make_model_using_centres():   _run('pex_make_model_using_centres')
def pex_modify_multiple():            _run('pex_modify_multiple')
def pex_peet_to_dynamo():             _run('pex_peet_to_dynamo')
def pex_peet_to_pytom():              _run('pex_peet_to_pytom')
def pex_plot_resolution():            _run('pex_plot_resolution')
def pex_plotback():                   _run('pex_plotback')
def pex_prm_prep():                   _run('pex_prm_prep')
def pex_randomise_pcles():            _run('pex_randomise_pcles')
def pex_remove_bad_symm_pcles():      _run('pex_remove_bad_symm_pcles')
def pex_remove_duplicates():          _run('pex_remove_duplicates')
def pex_remove_edge_pcles():          _run('pex_remove_edge_pcles')
def pex_remove_tube_duplicates():     _run('pex_remove_tube_duplicates')
def pex_revoke_when_complete():       _run('pex_revoke_when_complete')
def pex_run_sumcorr_for_tomo_SerEM_dev(): _run('pex_run_sumcorr_for_tomo_SerEM_dev')
def pex_scale_map():                  _run('pex_scale_map')
def pex_spherical_clean():            _run('pex_spherical_clean')
def pex_spherical_pick():             _run('pex_spherical_pick')
def pex_split_by_class():             _run('pex_split_by_class')
def pex_split_for_FSC():              _run('pex_split_for_FSC')
def pex_symm_checker():               _run('pex_symm_checker')
def pex_symmetrise():                 _run('pex_symmetrise')
def rex_make_subtomo_dose_file():     _run('rex_make_subtomo_dose_file')
def rex_star_to_newstack():           _run('rex_star_to_newstack')
def test_pex_unblur():                _run('test_pex_unblur')
