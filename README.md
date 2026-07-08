# peet_extensions (pex)

Written by Daven Vasishtan and Vojtech Prazak.
Extends the IMOD/PEET subvolume averaging pipeline.
Built on top of TEMPy (bundled in `TEMPy_py3/`).

---

## Requirements

- Python >= 3.8
- IMOD/PEET installed and on PATH (pex works alongside it, not as a replacement)
- conda or pip environment

Python dependencies are installed automatically (see [Installation](#installation)):

| Package | Purpose |
|---|---|
| numpy, scipy, biopython | TEMPy core |
| mrcfile | MRC/CCP4 map file I/O |
| matplotlib | plotting |
| scikit-image | image processing |
| scikit-learn | clustering |
| lxml | XML parsing |
| transforms3d, transformations | rotation/Euler angle utilities |
| PyCifRW | CIF/STAR file parsing (provides `CifFile`) |

---

## Installation

From the repo root, run:

```bash
pip install -e .
```

This installs all Python dependencies and makes all `pex_*` commands available on your PATH. The `-e` (editable) flag means changes to the source take effect immediately without reinstalling.

> **Note:** if you see a `SyntaxError` about TEMPy pointing to a file inside `anaconda/site-packages`, a stale manually-installed copy of TEMPy is shadowing the bundled version. Remove it with:
> ```bash
> rm -rf $(python -c "import TEMPy; import os; print(os.path.dirname(TEMPy.__file__))")
> ```
> Then reinstall: `pip install -e .`

---

## Running commands

All commands support `--help`:

```bash
pex_get_pentons --help
```

Commands can also be called directly by path (no PATH setup needed):

```bash
/path/to/repo/TEMPy_extensions_py3/pex_scripts/pex_get_pentons --help
```

---

## Command reference

### Particle picking

| Command | Description |
|---|---|
| `pex_get_pentons` | Pick pentons from an icosahedral capsid average |
| `pex_get_hexons` | Pick hexons from an icosahedral capsid average |
| `pex_get_icos_symm_pcles` | Pick particles using icosahedral symmetry |
| `pex_get_tetra` | Pick particles using tetrahedral symmetry |
| `pex_spherical_pick` | Place oriented particles on spheres from a model |
| `pex_linear_pick` | Place oriented particles along lines from a model |

### Particle cleaning

| Command | Description |
|---|---|
| `pex_clean_by_tilt_ang` | Remove/reset particles by tilt angle change |
| `pex_clean_by_planes` | Remove particles above/below planes in a model |
| `pex_clean_using_chim_markers` | Remove particles using a Chimera marker file |
| `pex_clean_using_imod_model` | Remove particles using an IMOD model file |
| `pex_remove_duplicates` | Remove overlapping particles (keeps highest CCC) |
| `pex_remove_tube_duplicates` | Remove duplicates along helices |
| `pex_remove_edge_pcles` | Remove particles too close to tomogram edges |
| `pex_remove_bad_symm_pcles` | Remove particles with symmetry partner issues |
| `pex_spherical_clean` | Remove particles outside a spherical radius |

### Particle manipulation

| Command | Description |
|---|---|
| `pex_modify_multiple` | Apply translation/rotation to particles |
| `pex_randomise_pcles` | Randomise particle angles up to a maximum |
| `pex_split_by_class` | Split particles by class assignment |
| `pex_split_for_FSC` | Split into even/odd halves for FSC |
| `pex_symmetrise` | Create new motive lists/models for symmetrisation |
| `pex_bin_model_files` | Bin particle positions for lower-resolution data |
| `pex_combine_model_files` | Combine model files per tomogram |
| `pex_combine_prm_files` | Combine multiple PEET parameter files into one |

### Visualisation and analysis

| Command | Description |
|---|---|
| `pex_make_chim_markers` | Make Chimera marker file from motl + model |
| `pex_make_chim_markers_from_prm` | Make Chimera marker file from a prm file |
| `pex_make_bilds` | Make `.bild` angular distribution file for Chimera |
| `pex_make_link_markers` | Make Chimera markers showing particle links |
| `pex_angular_distribution` | Get angle distribution relative to Euclidian planes |
| `pex_get_cccs` | Get cross-correlation values for an iteration |
| `pex_symm_checker` | Check C-symmetry of a map by rotational correlation |
| `pex_plot_resolution` | Plot FSC resolution curve |
| `pex_plotback` | Plot average back into tomogram |

### Map operations

| Command | Description |
|---|---|
| `pex_add_maps` | Add two MRC files together |
| `pex_scale_map` | Resample a map to a given voxel size |
| `pex_make_cylinder` | Create a solid cylinder MRC file |
| `pex_make_model_using_centres` | Make a model from a set of centres |

### Workflow utilities

| Command | Description |
|---|---|
| `pex_prm_prep` | Format initial PEET parameter file |
| `pex_revoke_when_complete` | Terminate Slurm allocation after PEET finishes |
| `pex_get_zeros_of_ctf` | Get CTF zero positions for given parameters |
| `pex_peet_to_dynamo` | Convert PEET output to Dynamo table format |
| `pex_peet_to_pytom` | Convert PEET output to PyTom XML format |

### RELION / format conversion

| Command | Description |
|---|---|
| `rex_star_to_newstack` | Convert Relion STAR file to newstack fileinlist |
| `rex_make_subtomo_dose_file` | Create cumulative dose file for Relion subtomo |
| `etomo_from_aretomo2` | Generate IMOD directives from header and run etomo |
| `make_fid_out_of_3d_model` | Create fiducial file from a 3D model |

---

## Known issues

Some scripts use deprecated scipy APIs (`scipy.ndimage.filters`, `scipy.ndimage.interpolation`) removed in scipy >= 1.12. If you see an `ImportError` from scipy, pin to a compatible version:

```bash
pip install "scipy>=1.7,<1.12"
```
