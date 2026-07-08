"""
End-to-end tests for pex commands that call library functions directly.
All tests that need real data use the data_dir / peet_dir fixtures from
conftest.py and are automatically skipped if the data is not present.
"""
import os
import pytest
import numpy as np
from xml.etree import ElementTree
from PEETModelParser import PEETmodel
from PEETMotiveList import PEETMotiveList


# ---------------------------------------------------------------------------
# pex_linear_pick
# ---------------------------------------------------------------------------

def test_linear_pick_produces_output_files(data_dir, tmp_path):
    from PEETPicker import linear_sampling
    m = PEETmodel(str(data_dir / 'two_points.mod'))
    out = str(tmp_path / 'lp_out')
    linear_sampling(m, 1.0, out, rand_orient=False)
    assert (tmp_path / 'lp_out.csv').exists()
    assert (tmp_path / 'lp_out.mod').exists()


def test_linear_pick_csv_has_particles(data_dir, tmp_path):
    from PEETPicker import linear_sampling
    m = PEETmodel(str(data_dir / 'two_points.mod'))
    out = str(tmp_path / 'lp_out')
    linear_sampling(m, 1.0, out, rand_orient=False)
    csv = PEETMotiveList(str(tmp_path / 'lp_out.csv'))
    assert len(csv) > 0


def test_linear_pick_is_deterministic(data_dir, tmp_path):
    """Non-random orientation mode must produce the same particle count twice."""
    from PEETPicker import linear_sampling
    m = PEETmodel(str(data_dir / 'two_points.mod'))

    out1 = str(tmp_path / 'run1_out')
    out2 = str(tmp_path / 'run2_out')
    linear_sampling(m, 1.0, out1, rand_orient=False)
    linear_sampling(m, 1.0, out2, rand_orient=False)

    n1 = len(PEETMotiveList(out1 + '.csv'))
    n2 = len(PEETMotiveList(out2 + '.csv'))
    assert n1 == n2


def test_linear_pick_creates_valid_cmm(data_dir, tmp_path):
    """pex_linear_pick internally calls csvfile_to_chim_markers; verify .cmm is valid XML."""
    from PEETPicker import linear_sampling
    from PEETParticleAnalysis import csvfile_to_chim_markers

    m = PEETmodel(str(data_dir / 'two_points.mod'))
    out = str(tmp_path / 'lp_out')
    linear_sampling(m, 1.0, out, rand_orient=False)

    cmm = out + '.cmm'
    csvfile_to_chim_markers(out + '.csv', out + '.mod', cmm)

    assert (tmp_path / 'lp_out.cmm').exists()
    tree = ElementTree.parse(cmm)
    markers = tree.findall('.//marker')
    assert len(markers) > 0


# ---------------------------------------------------------------------------
# pex_make_chim_markers — dedicated test using pre-existing lp output files
# ---------------------------------------------------------------------------

def test_make_chim_markers_from_existing_lp_files(data_dir, tmp_path):
    """pex_make_chim_markers: csv + mod → cmm with correct marker count."""
    from PEETParticleAnalysis import csvfile_to_chim_markers

    csv_path = str(data_dir / 'two_points_lp.csv')
    mod_path = str(data_dir / 'two_points_lp.mod')
    out_cmm = str(tmp_path / 'out.cmm')

    csvfile_to_chim_markers(csv_path, mod_path, out_cmm)

    assert (tmp_path / 'out.cmm').exists()
    tree = ElementTree.parse(out_cmm)
    markers = tree.findall('.//marker')
    n_particles = len(PEETMotiveList(csv_path))
    # Each particle produces a head + tail marker
    assert len(markers) == 2 * n_particles


# ---------------------------------------------------------------------------
# pex_spherical_pick
# ---------------------------------------------------------------------------

def test_spherical_pick_n1_particle_count(data_dir, tmp_path):
    """N=1 icosphere on a 2-point sphere model → 42 particles."""
    from PEETPicker import get_spherical_pick_from_model
    from PEETParticleAnalysis import csvfile_to_chim_markers

    m = PEETmodel(str(data_dir / 'two_points.mod'))
    out = str(tmp_path / 'sph_out')
    csv, mod = get_spherical_pick_from_model(m, 1, outfile=out)

    assert len(csv) == 42
    assert len(mod) == 42


def test_spherical_pick_writes_csv_and_mod(data_dir, tmp_path):
    from PEETPicker import get_spherical_pick_from_model

    m = PEETmodel(str(data_dir / 'two_points.mod'))
    out = str(tmp_path / 'sph_out')
    get_spherical_pick_from_model(m, 1, outfile=out)

    assert (tmp_path / 'sph_out.csv').exists()
    assert (tmp_path / 'sph_out.mod').exists()


def test_spherical_pick_cmm_is_valid_xml(data_dir, tmp_path):
    """Marker file produced for spherical pick is valid XML with markers."""
    from PEETPicker import get_spherical_pick_from_model
    from PEETParticleAnalysis import csvfile_to_chim_markers

    m = PEETmodel(str(data_dir / 'two_points.mod'))
    out = str(tmp_path / 'sph_out')
    get_spherical_pick_from_model(m, 1, outfile=out)
    csvfile_to_chim_markers(out + '.csv', out + '.mod', out + '.cmm')

    tree = ElementTree.parse(out + '.cmm')
    assert len(tree.findall('.//marker')) == 2 * 42


# ---------------------------------------------------------------------------
# pex_bin_model_files
# ---------------------------------------------------------------------------

def test_bin_model_preserves_particle_count(data_dir, tmp_path, cwd_run1):
    """Binning by 2 should not change particle count."""
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    prm.bin_all(2.0, 1, str(tmp_path / 'bin2'), writeprm=False)

    out_csvs = list((tmp_path / 'bin2').glob('*.csv'))
    assert len(out_csvs) == 1
    binned = PEETMotiveList(str(out_csvs[0]))
    assert len(binned) == 258


def test_bin_model_halves_coordinates(data_dir, peet_dir, tmp_path, cwd_run1):
    """Binning by 2 should halve all particle positions."""
    from PEETPRMParser import PEETPRMFile

    # The mod file referenced in the prm (../r05_1_*.mod) lives in peet_dir
    orig_mod = PEETmodel(str(peet_dir / 'r05_1_MOTL_Tom1_Iter7_remdup_2.0.mod'))
    orig_pts = orig_mod.get_all_points().copy()

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    prm.bin_all(2.0, 1, str(tmp_path / 'bin2'), writeprm=False, offset_mv=False)

    out_mods = list((tmp_path / 'bin2').glob('*.mod'))
    binned_mod = PEETmodel(str(out_mods[0]))
    binned_pts = binned_mod.get_all_points()

    np.testing.assert_allclose(binned_pts, orig_pts / 2.0, atol=1e-3)


# ---------------------------------------------------------------------------
# pex_randomise_pcles
# ---------------------------------------------------------------------------

def test_randomise_preserves_particle_count(data_dir, tmp_path, cwd_run1):
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    prm.randomise_pcles(1, 5.0, str(tmp_path / 'rand'), writeprm=False)

    out_csvs = list((tmp_path / 'rand').glob('*.csv'))
    assert len(out_csvs) == 1
    rand = PEETMotiveList(str(out_csvs[0]))
    assert len(rand) == 258


def test_randomise_changes_angles(data_dir, tmp_path, cwd_run1):
    """Randomisation should change at least some particle angles."""
    from PEETPRMParser import PEETPRMFile

    orig = PEETMotiveList(str(data_dir / 'run1' / 'r_MOTL_Tom1_Iter1.csv'))
    orig_angles = np.array([orig.get_angles_by_list_index(i) for i in range(len(orig))])

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    prm.randomise_pcles(1, 180.0, str(tmp_path / 'rand'), writeprm=False)

    out_csvs = list((tmp_path / 'rand').glob('*.csv'))
    rand = PEETMotiveList(str(out_csvs[0]))
    rand_angles = np.array([rand.get_angles_by_list_index(i) for i in range(len(rand))])

    assert not np.allclose(orig_angles, rand_angles, atol=1e-6)


# ---------------------------------------------------------------------------
# pex_split_for_FSC
# ---------------------------------------------------------------------------

def test_split_for_fsc_counts_sum_to_total(data_dir, tmp_path, cwd_run1):
    """Even + odd halves must together account for all 258 particles."""
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    # split_by_classID uses os.mkdir (not makedirs), so pass flat dirs under tmp_path
    cls1_dir = str(tmp_path / 'fsc_cls1')
    cls2_dir = str(tmp_path / 'fsc_cls2')
    prm.split_by_classID(1, cls1_dir, classes=[1], splitForFSC=True, writeprm=False)
    prm.split_by_classID(1, cls2_dir, classes=[2], splitForFSC=True, writeprm=False)

    cls1_csvs = list((tmp_path / 'fsc_cls1').glob('*.csv'))
    cls2_csvs = list((tmp_path / 'fsc_cls2').glob('*.csv'))
    assert len(cls1_csvs) == 1 and len(cls2_csvs) == 1

    n1 = len(PEETMotiveList(str(cls1_csvs[0])))
    n2 = len(PEETMotiveList(str(cls2_csvs[0])))
    assert n1 + n2 == 258


def test_split_for_fsc_halves_are_roughly_equal(data_dir, tmp_path, cwd_run1):
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    cls1_dir = str(tmp_path / 'fsc_cls1')
    cls2_dir = str(tmp_path / 'fsc_cls2')
    prm.split_by_classID(1, cls1_dir, classes=[1], splitForFSC=True, writeprm=False)
    prm.split_by_classID(1, cls2_dir, classes=[2], splitForFSC=True, writeprm=False)

    n1 = len(PEETMotiveList(str(list((tmp_path / 'fsc_cls1').glob('*.csv'))[0])))
    n2 = len(PEETMotiveList(str(list((tmp_path / 'fsc_cls2').glob('*.csv'))[0])))
    # Allow ±1 for odd totals
    assert abs(n1 - n2) <= 1


# ---------------------------------------------------------------------------
# pex_split_by_class
# ---------------------------------------------------------------------------

def test_split_by_class_extracts_all_when_single_class(data_dir, tmp_path, cwd_run1):
    """iter1 has class=1 for all 258 particles; extracting class 1 keeps them all."""
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    prm.split_by_classID(1, str(tmp_path / 'cls'), classes=[1], splitForFSC=False, writeprm=False)

    out_csvs = list((tmp_path / 'cls').glob('*.csv'))
    assert len(out_csvs) == 1
    assert len(PEETMotiveList(str(out_csvs[0]))) == 258


def test_split_by_class_absent_class_produces_no_output(data_dir, tmp_path, cwd_run1):
    """Extracting a non-existent class ID should produce no output files."""
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    prm.split_by_classID(1, str(tmp_path / 'cls99'), classes=[99], splitForFSC=False, writeprm=False)

    out_csvs = list((tmp_path / 'cls99').glob('*.csv'))
    assert len(out_csvs) == 0


# ---------------------------------------------------------------------------
# pex_angular_distribution
# ---------------------------------------------------------------------------

def test_angular_distribution_shape(data_dir, cwd_run1):
    """get_angle_distributions returns a (n_particles, 3) array."""
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    dists = prm.get_angle_distributions(1)

    assert dists.ndim == 2
    assert dists.shape == (258, 3)


def test_angular_distribution_writes_file(data_dir, tmp_path, cwd_run1):
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    out = str(tmp_path / 'ang_dist.txt')
    prm.get_angle_distributions(1, outfile=out)

    assert (tmp_path / 'ang_dist.txt').exists()
    data = np.loadtxt(out)
    assert data.shape == (258, 3)


# ---------------------------------------------------------------------------
# pex_remove_duplicates (via direct function — already tested)
# ---------------------------------------------------------------------------

def test_remove_duplicates_preserves_all_at_zero_distance(peet_dir):
    """At max_dist=0, no particles should be considered overlapping."""
    from PEETParticleCleanup import remove_duplicates

    csv_path = str(peet_dir / 'r05_1_MOTL_Tom1_Iter7_remdup_2.0.csv')
    mod_path = str(peet_dir / 'r05_1_MOTL_Tom1_Iter7_remdup_2.0.mod')

    out_csv, out_mod, n_before, n_after = remove_duplicates(csv_path, mod_path, max_dist=0.0, verbose=False)
    assert n_before == 258
    assert n_after == 258
    assert len(out_csv) == 258
    assert len(out_mod) == 258


def test_remove_duplicates_output_matches_input_length(peet_dir):
    """Sanity check: output csv and mod lengths match each other."""
    from PEETParticleCleanup import remove_duplicates

    csv_path = str(peet_dir / 'r05_1_MOTL_Tom1_Iter7_remdup_2.0.csv')
    mod_path = str(peet_dir / 'r05_1_MOTL_Tom1_Iter7_remdup_2.0.mod')

    out_csv, out_mod, n_before, n_after = remove_duplicates(csv_path, mod_path, max_dist=0.0, verbose=False)
    assert len(out_csv) == len(out_mod)
    assert n_before == n_after


# ---------------------------------------------------------------------------
# pex_remove_edge_pcles
# ---------------------------------------------------------------------------

def test_remove_edge_pcles_zero_dist_keeps_all(data_dir, tmp_path, cwd_run1):
    """edge_dist=0 keeps every particle (none can be closer than 0 to an edge)."""
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    result = prm.remove_offedge_pcles(1, 0.0, str(tmp_path / 'edge0'), writeprm=False, verbose=False)

    out_csvs = list((tmp_path / 'edge0').glob('*.csv'))
    assert len(out_csvs) == 1
    assert len(PEETMotiveList(str(out_csvs[0]))) == 258


def test_remove_edge_pcles_large_dist_removes_all(data_dir, tmp_path, cwd_run1):
    """edge_dist larger than any tomogram dimension removes all particles."""
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    # 10000 px is larger than any real tomogram
    prm.remove_offedge_pcles(1, 10000.0, str(tmp_path / 'edge_big'), writeprm=False, verbose=False)

    # When all particles are removed the function prints a message and creates no files
    out_csvs = list((tmp_path / 'edge_big').glob('*.csv'))
    assert len(out_csvs) == 0


# ---------------------------------------------------------------------------
# pex_make_chim_markers_from_prm — prm → per-motl cmm
# ---------------------------------------------------------------------------

def test_make_chim_markers_from_prm_writes_cmm(data_dir, peet_dir, tmp_path):
    """make_chim_markers on a synthetic prm: output cmm lands in tmp_path, marker count == 2×258."""
    import shutil
    from PEETPRMParser import PEETPRMFile

    # Copy motl so the derived output path (motl[:-4]+'_markers.cmm') lands in tmp_path
    motl_copy = str(tmp_path / 'r_MOTL_Tom1_Iter1.csv')
    shutil.copy(str(peet_dir / 'run1' / 'r_MOTL_Tom1_Iter1.csv'), motl_copy)

    prm = PEETPRMFile('')
    prm.prm_dict = {
        'fnOutput': 'r',
        'initMOTL': [motl_copy],
        'fnModParticle': [str(peet_dir / 'r05_1_MOTL_Tom1_Iter7_remdup_2.0.mod')],
    }
    prm.make_chim_markers(0)

    cmm = tmp_path / 'r_MOTL_Tom1_Iter1_markers.cmm'
    assert cmm.exists()
    tree = ElementTree.parse(str(cmm))
    markers = tree.findall('.//marker')
    assert len(markers) == 2 * 258


# ---------------------------------------------------------------------------
# pex_make_bilds — prm + iter → bild arrow file
# ---------------------------------------------------------------------------

def test_make_bilds_writes_bild(data_dir, tmp_path, cwd_run1):
    from bildmaker import make_bild_from_prm

    outfile = str(tmp_path / 'test.bild')
    make_bild_from_prm(str(data_dir / 'run1' / 'r.prm'), 1, outfile, combined=True)

    assert (tmp_path / 'test.bild').exists()
    with open(outfile) as f:
        first_line = f.readline()
    assert first_line.startswith('.color')


# ---------------------------------------------------------------------------
# pex_get_cccs — prm + iter → per-particle CCC values
# ---------------------------------------------------------------------------

def test_get_cccs_returns_correct_count(data_dir, cwd_run1):
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    cccs = prm.get_all_ccc(1)

    assert len(cccs) == 258
    assert np.all(np.isfinite(cccs))


def test_get_cccs_writes_file(data_dir, tmp_path, cwd_run1):
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    outfile = str(tmp_path / 'cccs.txt')
    prm.get_all_ccc(1, outfile=outfile)

    assert (tmp_path / 'cccs.txt').exists()
    data = np.loadtxt(outfile)
    assert len(data) == 258


# ---------------------------------------------------------------------------
# pex_get_zeros_of_ctf — pure CTF math, no data needed
# ---------------------------------------------------------------------------

def test_get_zeros_of_ctf_returns_valid_curve():
    from ctf_calc import ctf

    r, curve = ctf(1.5, 300, 2.7, 9.1)

    assert r.shape == curve.shape
    assert len(r) > 0
    assert np.all(np.isfinite(curve))
    assert np.any(curve != 0)


# ---------------------------------------------------------------------------
# pex_make_cylinder — pure shape generation, no data needed
# ---------------------------------------------------------------------------

def test_make_cylinder_has_correct_box_size():
    from make_shapes import make_cylinder

    cyl = make_cylinder([50, 50, 50], 20, 40)

    assert cyl.fullMap.shape == (50, 50, 50)
    assert cyl.fullMap.sum() > 0


# ---------------------------------------------------------------------------
# pex_make_link_markers — csv + mod → cmm with distance links
# ---------------------------------------------------------------------------

def test_make_link_markers_writes_cmm(data_dir, tmp_path):
    from PEETParticleAnalysis import pcle_links

    motl = PEETMotiveList(str(data_dir / 'two_points_lp.csv'))
    mod = PEETmodel(str(data_dir / 'two_points_lp.mod'))
    outfile = str(tmp_path / 'links.cmm')

    pcle_links(motl, mod, max_dist=200.0, apix=9.105, outfile=outfile)

    assert (tmp_path / 'links.cmm').exists()
    tree = ElementTree.parse(outfile)
    markers = tree.findall('.//marker')
    assert len(markers) > 0


# ---------------------------------------------------------------------------
# pex_symmetrise — C2-symmetrise a motive list
# ---------------------------------------------------------------------------

def test_symmetrise_c2_creates_two_csvs(data_dir, tmp_path, cwd_run1):
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    prm.symmetrise('C2', 1, str(tmp_path / 'symm'), writeprm=False)

    out_csvs = sorted((tmp_path / 'symm').glob('*.csv'))
    assert len(out_csvs) == 2
    for csv_path in out_csvs:
        assert len(PEETMotiveList(str(csv_path))) == 258


def test_symmetrise_c2_rotates_angles(data_dir, tmp_path, cwd_run1):
    """The 0° and 180° symmetric copies must differ in at least some angles."""
    from PEETPRMParser import PEETPRMFile

    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    prm.symmetrise('C2', 1, str(tmp_path / 'symm'), writeprm=False)

    csvs = sorted((tmp_path / 'symm').glob('*.csv'))
    motl_0 = PEETMotiveList(str(csvs[0]))
    motl_180 = PEETMotiveList(str(csvs[1]))

    angles_0 = np.array([motl_0.get_angles_by_list_index(i) for i in range(10)])
    angles_180 = np.array([motl_180.get_angles_by_list_index(i) for i in range(10)])
    assert not np.allclose(angles_0, angles_180, atol=1e-3)


# ---------------------------------------------------------------------------
# pex_spherical_clean — filter particles outside sphere radius (synthetic)
# ---------------------------------------------------------------------------

def test_spherical_clean_large_tolerance_keeps_all(data_dir, tmp_path):
    """rad_extra much larger than sphere radius → all 42 particles kept."""
    from PEETPicker import get_spherical_pick_from_model
    from PEETSphericalClean import clean_spheres

    m = PEETmodel(str(data_dir / 'two_points.mod'))
    out = str(tmp_path / 'sph')
    get_spherical_pick_from_model(m, 1, outfile=out)

    _, _, in_num, out_num = clean_spheres(out + '.csv', out + '.mod',
                                           fixed_rad=None, rad_extra=10000, verbose=False)
    assert in_num == 42
    assert out_num == 0


def test_spherical_clean_tiny_radius_removes_all(data_dir, tmp_path):
    """fixed_rad=1, rad_extra=0 → no particle can lie on the 1-px shell → all removed."""
    from PEETPicker import get_spherical_pick_from_model
    from PEETSphericalClean import clean_spheres

    m = PEETmodel(str(data_dir / 'two_points.mod'))
    out = str(tmp_path / 'sph')
    get_spherical_pick_from_model(m, 1, outfile=out)

    _, _, in_num, out_num = clean_spheres(out + '.csv', out + '.mod',
                                           fixed_rad=1, rad_extra=0, verbose=False)
    assert in_num == 0
    assert out_num == 42


# ---------------------------------------------------------------------------
# pex_modify_multiple — translate/rotate all particles (uses IMOD modifyMotiveList)
# ---------------------------------------------------------------------------

def test_modify_multiple_zero_translation_preserves_count(data_dir, tmp_path, cwd_run1):
    """Zero translation and zero rotation must leave particle count unchanged."""
    from modify_multiple_lib import modify_prm

    outdir = str(tmp_path / 'mod_out')
    modify_prm(str(data_dir / 'run1' / 'r.prm'), 1,
               [0, 0, 0], [0, 0, 0], outdir, max_cores=1)

    out_csvs = list((tmp_path / 'mod_out').glob('*.csv'))
    assert len(out_csvs) == 1
    assert len(PEETMotiveList(str(out_csvs[0]))) == 258


# ---------------------------------------------------------------------------
# pex_remove_tube_duplicates — expand then deduplicate (tests library calls)
# ---------------------------------------------------------------------------

def test_remove_tube_duplicates_expand_step(data_dir, tmp_path, cwd_run1):
    """Expand step (modify_prm with width=[10,0,0]) must preserve particle count."""
    from modify_multiple_lib import modify_prm

    expand_dir = str(tmp_path / 'expand')
    modify_prm(str(data_dir / 'run1' / 'r.prm'), 1,
               [10, 0, 0], [0, 0, 0], expand_dir, max_cores=1)

    out_csvs = list((tmp_path / 'expand').glob('*.csv'))
    assert len(out_csvs) == 1
    assert len(PEETMotiveList(str(out_csvs[0]))) == 258


def test_remove_tube_duplicates_dedup_step(data_dir, tmp_path, cwd_run1):
    """clean_pcles at max_dist=0 (no particles within 0 Å of each other) keeps all 258."""
    from modify_multiple_lib import modify_prm
    from PEETPRMParser import PEETPRMFile
    import glob

    expand_dir = str(tmp_path / 'expand')
    modify_prm(str(data_dir / 'run1' / 'r.prm'), 1,
               [10, 0, 0], [0, 0, 0], expand_dir, max_cores=1)

    expanded_prm_path = glob.glob(expand_dir + '/*.prm')[0]
    expanded_prm = PEETPRMFile(expanded_prm_path)
    dedup_dir = str(tmp_path / 'dedup')
    expanded_prm.clean_pcles(0, 0, dedup_dir, writeprm=False, verbose=False)

    out_csvs = list((tmp_path / 'dedup').glob('*.csv'))
    assert len(out_csvs) == 1
    assert len(PEETMotiveList(str(out_csvs[0]))) == 258


# ---------------------------------------------------------------------------
# pex_peet_to_dynamo — write Dynamo table without particle extraction
# ---------------------------------------------------------------------------

def test_peet_to_dynamo_writes_table(data_dir, tmp_path, cwd_run1):
    """extpcles=False skips tomogram access; particles.tbl must have 258 rows."""
    from PEET2Dynamo import extract_pcles_for_dynamo_from_prmfile

    outdir = str(tmp_path / 'dynamo')
    os.makedirs(outdir)
    extract_pcles_for_dynamo_from_prmfile(
        str(data_dir / 'run1' / 'r.prm'), 1, [32, 32, 32], outdir, extpcles=False
    )

    tbl = tmp_path / 'dynamo' / 'particles.tbl'
    assert tbl.exists()
    data = np.loadtxt(str(tbl))
    assert data.shape[0] == 258
    assert np.all(np.isfinite(data))


# ---------------------------------------------------------------------------
# pex_peet_to_pytom — write PyTom XML lists without particle extraction
# ---------------------------------------------------------------------------

def test_peet_to_pytom_writes_xml(data_dir, cwd_run1):
    """extpcles=False writes per-tom and combined XML; files must parse and contain particles."""
    from PEET2Pytom import extract_pcles_using_PEET_PRM

    output_files = ['tom001_pcle_list.xml', 'all_pcle_lists.xml', 'job.xml']
    try:
        extract_pcles_using_PEET_PRM(
            str(data_dir / 'run1' / 'r.prm'), 1, [32, 32, 32], extpcles=False
        )
        assert os.path.exists('tom001_pcle_list.xml')
        assert os.path.exists('all_pcle_lists.xml')
        assert os.path.exists('job.xml')

        tree = ElementTree.parse('tom001_pcle_list.xml')
        pcles = tree.findall('.//Particle')
        assert len(pcles) == 258
    finally:
        for fname in output_files:
            if os.path.exists(fname):
                os.remove(fname)


# ---------------------------------------------------------------------------
# pex_clean_by_planes
# ---------------------------------------------------------------------------

def test_clean_by_planes_keeps_all_when_plane_is_above(data_dir, tmp_path):
    """Plane at z=99999 with remove='above': all particles are below the plane → all kept."""
    from clean_by_plane import clean_by_y_orth_plane

    csv = PEETMotiveList(str(data_dir / 'two_points_lp.csv'))
    mod = PEETmodel(str(data_dir / 'two_points_lp.mod'))
    n_total = len(csv)

    plane_mod = PEETmodel()
    plane_mod.add_point(0, 0, [0.0, 0.0, 99999.0])
    plane_mod.add_point(0, 0, [100.0, 0.0, 99999.0])
    plane_modfile = str(tmp_path / 'plane.mod')
    plane_mod.write_model(plane_modfile)

    new_csv, new_mod = clean_by_y_orth_plane(csv, mod, plane_modfile, remove='above')
    assert len(new_csv) == n_total
    assert len(new_mod) == n_total


def test_clean_by_planes_removes_all_when_plane_is_above(data_dir, tmp_path):
    """Plane at z=99999 with remove='below': z_lim always > particle z → none kept."""
    from clean_by_plane import clean_by_y_orth_plane

    csv = PEETMotiveList(str(data_dir / 'two_points_lp.csv'))
    mod = PEETmodel(str(data_dir / 'two_points_lp.mod'))

    plane_mod = PEETmodel()
    plane_mod.add_point(0, 0, [0.0, 0.0, 99999.0])
    plane_mod.add_point(0, 0, [100.0, 0.0, 99999.0])
    plane_modfile = str(tmp_path / 'plane.mod')
    plane_mod.write_model(plane_modfile)

    new_csv, new_mod = clean_by_y_orth_plane(csv, mod, plane_modfile, remove='below')
    assert len(new_csv) == 0
    assert len(new_mod) == 0


# ---------------------------------------------------------------------------
# pex_combine_model_files
# ---------------------------------------------------------------------------

def test_combine_model_files_single_tomo_preserves_count(data_dir, tmp_path, cwd_run1):
    """With 1 tomogram, combined output equals the single input motive list."""
    from PEETPRMParser import PEETPRMFile
    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    outdir = str(tmp_path / 'combined')
    prm.combine_model_files(1, outdir, writeprm=False)
    out_csvs = list((tmp_path / 'combined').glob('*.csv'))
    assert len(out_csvs) == 1
    assert len(PEETMotiveList(str(out_csvs[0]))) == 258


# ---------------------------------------------------------------------------
# pex_combine_prm_files
# ---------------------------------------------------------------------------

def test_combine_prm_files_doubles_tomogram_count(data_dir, tmp_path, cwd_run1):
    """Combining a prm with itself produces a new prm with 2 tomogram entries."""
    from PEETPRMParser import PEETPRMFile
    p = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    combined = p.deepcopy()
    for key in ('initMOTL', 'fnModParticle', 'fnVolume', 'tiltRange'):
        combined.prm_dict[key] = p.prm_dict[key] + p.prm_dict[key]
    out = str(tmp_path / 'combined.prm')
    combined.write_prm_file(out)
    reread = PEETPRMFile(out)
    assert len(reread.prm_dict['fnVolume']) == 2
    assert len(reread.prm_dict['initMOTL']) == 2


# ---------------------------------------------------------------------------
# pex_get_pentons
# ---------------------------------------------------------------------------

def test_get_pentons_creates_12x_particles(data_dir, tmp_path, cwd_run1):
    """Each input particle → 12 pentons (icosahedral vertices)."""
    from PEETPRMParser import PEETPRMFile
    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    outdir = str(tmp_path / 'pentons')
    prm.get_pentons(1, 100.0, outdir, writeprm=False)
    csvs = list((tmp_path / 'pentons').glob('*.csv'))
    assert len(csvs) == 1
    n_out = len(PEETMotiveList(str(csvs[0])))
    assert n_out == 12 * 258
    mods = list((tmp_path / 'pentons').glob('*.mod'))
    assert len(PEETmodel(str(mods[0])).get_all_points()) == n_out


# ---------------------------------------------------------------------------
# pex_get_hexons
# ---------------------------------------------------------------------------

def test_get_hexons_creates_multiple_particles(data_dir, tmp_path, cwd_run1):
    """Each input particle → several hexons; output is a multiple of input count."""
    from PEETPRMParser import PEETPRMFile
    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    outdir = str(tmp_path / 'hexons')
    prm.get_hexons(1, 100.0, 2, outdir, writeprm=False)
    csvs = list((tmp_path / 'hexons').glob('*.csv'))
    assert len(csvs) == 1
    n_out = len(PEETMotiveList(str(csvs[0])))
    assert n_out > 258
    assert n_out % 258 == 0
    mods = list((tmp_path / 'hexons').glob('*.mod'))
    assert len(PEETmodel(str(mods[0])).get_all_points()) == n_out


# ---------------------------------------------------------------------------
# pex_get_tetra
# ---------------------------------------------------------------------------

def test_get_tetra_creates_4x_particles(data_dir, tmp_path, cwd_run1):
    """Each input particle → 4 tetrahedral vertices."""
    from PEETPRMParser import PEETPRMFile
    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    outdir = str(tmp_path / 'tetra')
    prm.get_tetra(1, 100.0, outdir, writeprm=False)
    csvs = list((tmp_path / 'tetra').glob('*.csv'))
    assert len(csvs) == 1
    n_out = len(PEETMotiveList(str(csvs[0])))
    assert n_out == 4 * 258
    mods = list((tmp_path / 'tetra').glob('*.mod'))
    assert len(PEETmodel(str(mods[0])).get_all_points()) == n_out


# ---------------------------------------------------------------------------
# pex_get_icos_symm_pcles
# ---------------------------------------------------------------------------

def test_get_icos_symm_pcles_creates_60x_particles(data_dir, tmp_path, cwd_run1):
    """Each input particle → 60 icosahedral symmetry-related positions."""
    from PEETPRMParser import PEETPRMFile
    from Vector import Vector
    prm = PEETPRMFile(str(data_dir / 'run1' / 'r.prm'))
    outdir = str(tmp_path / 'icos')
    prm.get_general_icos_pos(1, Vector.fromlist([1.0, 0.0, 0.0]), outdir, writeprm=False)
    csvs = list((tmp_path / 'icos').glob('*.csv'))
    assert len(csvs) == 1
    n_out = len(PEETMotiveList(str(csvs[0])))
    assert n_out == 60 * 258
    mods = list((tmp_path / 'icos').glob('*.mod'))
    assert len(PEETmodel(str(mods[0])).get_all_points()) == n_out


# ---------------------------------------------------------------------------
# pex_plot_resolution
# ---------------------------------------------------------------------------

def test_plot_resolution_fsc_returns_valid_curve(peet_dir):
    """fsc() on run2 half-maps returns finite FSC values in [-1, 1]."""
    import mrcfile
    from PEET_matlab_trans_fn import fsc
    ref1 = str(peet_dir / 'run2' / 'r2_Ref1.mrc')
    ref2 = str(peet_dir / 'run2' / 'r2_Ref2.mrc')
    with mrcfile.open(ref1, mode='r', permissive=True) as m:
        v1 = m.data.astype(np.float32).copy()
    with mrcfile.open(ref2, mode='r', permissive=True) as m:
        v2 = m.data.astype(np.float32).copy()
    scores, freqs = fsc(v1, v2)
    assert len(scores) > 0
    assert np.all(np.isfinite(scores))
    assert np.all(scores >= -1.0) and np.all(scores <= 1.0)


# ---------------------------------------------------------------------------
# pex_symm_checker
# ---------------------------------------------------------------------------

def test_symm_checker_runs_without_error(peet_dir):
    """symm_checker on the run1 average with step=30 (12 rotations) completes cleanly."""
    import matplotlib
    matplotlib.use('Agg')
    from symm_checker import symm_checker
    symm_checker(
        str(peet_dir / 'run1' / 'r_AvgVol_1P0258.mrc'),
        outfile='',
        verbose=False,
        step=30,
    )


# ---------------------------------------------------------------------------
# pex_add_maps
# ---------------------------------------------------------------------------

def test_add_maps_output_exists(tmp_path):
    import mrcfile
    from MapParser_f32_new import MapParser
    data = np.ones([20, 20, 20], dtype=np.float32)
    for name in ('a.mrc', 'b.mrc'):
        with mrcfile.new(str(tmp_path / name), overwrite=True) as mrc:
            mrc.set_data(data)
            mrc.voxel_size = 1.0
    m1 = MapParser.readMRC(str(tmp_path / 'a.mrc'))
    m2 = MapParser.readMRC(str(tmp_path / 'b.mrc'))
    m1.fullMap = m1.fullMap + m2.fullMap
    outfile = str(tmp_path / 'out.mrc')
    m1.write_to_MRC_file(outfile)
    assert (tmp_path / 'out.mrc').exists()


def test_add_maps_values_are_correct(tmp_path):
    """Adding maps with weights w1=1, w2=2 should give voxel values of 3.0."""
    import mrcfile
    from MapParser_f32_new import MapParser
    with mrcfile.new(str(tmp_path / 'a.mrc'), overwrite=True) as mrc:
        mrc.set_data(np.ones([20, 20, 20], dtype=np.float32))
        mrc.voxel_size = 1.0
    with mrcfile.new(str(tmp_path / 'b.mrc'), overwrite=True) as mrc:
        mrc.set_data(np.ones([20, 20, 20], dtype=np.float32) * 2.0)
        mrc.voxel_size = 1.0
    m1 = MapParser.readMRC(str(tmp_path / 'a.mrc'))
    m2 = MapParser.readMRC(str(tmp_path / 'b.mrc'))
    m1.fullMap = m1.fullMap * 1.0 + m2.fullMap * 1.0
    outfile = str(tmp_path / 'out.mrc')
    m1.write_to_MRC_file(outfile)
    with mrcfile.open(outfile, permissive=True) as mrc:
        assert np.allclose(mrc.data, 3.0, atol=1e-4)


# ---------------------------------------------------------------------------
# pex_scale_map
# ---------------------------------------------------------------------------

def test_scale_map_changes_voxel_size(tmp_path):
    import mrcfile
    from MapParser_f32_new import MapParser
    with mrcfile.new(str(tmp_path / 'in.mrc'), overwrite=True) as mrc:
        mrc.set_data(np.random.default_rng(0).random([20, 20, 20]).astype(np.float32))
        mrc.voxel_size = 1.0
    m = MapParser.readMRC(str(tmp_path / 'in.mrc'))
    out = m.resample_by_apix(2.0)
    outfile = str(tmp_path / 'scaled.mrc')
    out.write_to_MRC_file(outfile)
    with mrcfile.open(outfile, permissive=True) as mrc:
        assert abs(float(mrc.voxel_size.x) - 2.0) < 0.5


# ---------------------------------------------------------------------------
# pex_plotback
# ---------------------------------------------------------------------------

def test_plotback_writes_output_mrc(data_dir, tmp_path):
    """replace_pcles on two_points_lp particles (64) with synthetic tiny average."""
    import mrcfile
    from replace_pcles import replace_pcles

    # Synthetic 8×8×8 average with a non-zero centre voxel
    ave_path = str(tmp_path / 'ave.mrc')
    ave_data = np.zeros([8, 8, 8], dtype=np.float32)
    ave_data[4, 4, 4] = 1.0
    with mrcfile.new(ave_path, overwrite=True) as mrc:
        mrc.set_data(ave_data)
        mrc.voxel_size = 9.105

    outfile = str(tmp_path / 'plotback.mrc')
    # two_points_lp positions: x≈49-65, y≈743-804, z=102 → tomo [100, 810, 110]
    replace_pcles(
        ave_path,
        [100, 810, 110],
        str(data_dir / 'two_points_lp.csv'),
        str(data_dir / 'two_points_lp.mod'),
        outfile,
    )
    assert (tmp_path / 'plotback.mrc').exists()
    with mrcfile.open(outfile, permissive=True) as mrc:
        assert mrc.data.max() > 0


# ---------------------------------------------------------------------------
# Deferred / known-broken (placeholder)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# pex_make_model_using_centres
# ---------------------------------------------------------------------------

def test_make_model_using_centres_maps_each_pick_to_nearest_centre(data_dir):
    """Each pick particle is assigned a two-point contour (pick + nearest centre); no outfile → no stalkInit binary needed."""
    from make_model_using_centres_new import make_model
    pick_mod = str(data_dir / 'two_points_lp.mod')
    cen_mod  = str(data_dir / 'one_point.mod')
    result = make_model(pick_mod, cen_mod, outfile=None)
    n_pick = len(PEETmodel(pick_mod).get_all_points())
    assert len(result) == 2 * n_pick


# ---------------------------------------------------------------------------
# pex_clean_by_tilt_ang
# ---------------------------------------------------------------------------

def test_clean_by_tilt_ang_large_tolerance_keeps_all(data_dir, tmp_path):
    """Comparing a motive list to itself gives Δangle=0; max_ang=999 keeps everything."""
    from PEETParticleCleanup import clean_pcles_by_tilt_ang_change
    csv_f = str(data_dir / 'two_points_lp.csv')
    mod_f = str(data_dir / 'two_points_lp.mod')
    _, _, bef, aft = clean_pcles_by_tilt_ang_change(
        csv_f, mod_f, csv_f, max_ang=999, outfile=str(tmp_path / 'out'), verbose=False)
    assert aft == bef


def test_clean_by_tilt_ang_zero_tolerance_removes_all(data_dir, tmp_path):
    """max_ang=0 keeps condition 0 < 0 is False → removes all particles."""
    from PEETParticleCleanup import clean_pcles_by_tilt_ang_change
    csv_f = str(data_dir / 'two_points_lp.csv')
    mod_f = str(data_dir / 'two_points_lp.mod')
    _, _, _, aft = clean_pcles_by_tilt_ang_change(
        csv_f, mod_f, csv_f, max_ang=0, outfile='', verbose=False)
    assert aft == 0


# ---------------------------------------------------------------------------
# pex_remove_bad_symm_pcles (no-partners case)
# ---------------------------------------------------------------------------

def test_remove_bad_symm_pcles_no_partners_removes_all(data_dir):
    """t_tol=0.001 < min inter-particle dist (1px) → every particle has 0 C2 partners → all flagged bad."""
    from check_symm_pcles import check_symm_pcles
    _, _, _, bef, aft = check_symm_pcles(
        str(data_dir / 'two_points_lp.csv'),
        str(data_dir / 'two_points_lp.mod'),
        sym=2, t_tol=0.001, r_tol=999, outfile='', verbose=False)
    assert aft == 0


@pytest.mark.skip(reason="pex_remove_bad_symm_pcles positive case needs C-symmetric particle data")
def test_remove_bad_symm_pcles_placeholder():
    pass


# ---------------------------------------------------------------------------
# pex_clean_using_imod_model
# ---------------------------------------------------------------------------

@pytest.mark.skip(reason="pex_clean_using_imod_model: clean_prm_using_model is defined inside main() — not importable")
def test_clean_using_imod_model_placeholder():
    pass


# ---------------------------------------------------------------------------
# pex_clean_using_chim_markers
# ---------------------------------------------------------------------------

def test_clean_using_chim_markers_all_heads_present_keeps_all(data_dir, tmp_path):
    """Synthetic CMM with all particle IDs and b != '1' → all particles kept."""
    from PEETParticleCleanup import clean_using_marker_file
    n = len(PEETMotiveList(str(data_dir / 'two_points_lp.csv')))
    lines = ['<marker_set name="test">']
    for i in range(n):
        lines.append(f'  <marker id="{i}" x="0" y="0" z="0" r="1" g="0" b="0.5"/>')
    lines.append('</marker_set>')
    cmm = tmp_path / 'test.cmm'
    cmm.write_text('\n'.join(lines))
    _, _, bef, aft = clean_using_marker_file(
        str(data_dir / 'two_points_lp.csv'),
        str(data_dir / 'two_points_lp.mod'),
        str(cmm), outfile='', verbose=False)
    assert aft == n


@pytest.mark.skip(reason="rex_make_subtomo_dose_file needs IMOD extracttilts binary + tilt stack")
def test_make_subtomo_dose_file_placeholder():
    pass


@pytest.mark.skip(reason="rex_star_to_newstack needs CifFile module + Relion .star file")
def test_star_to_newstack_placeholder():
    pass


@pytest.mark.skip(reason="pex_revoke_when_complete is a Slurm HPC watcher — untestable without cluster")
def test_revoke_when_complete_placeholder():
    pass


@pytest.mark.skip(reason="pex_run_sumcorr_for_tomo_SerEM_dev needs MotionCor2 — skip all _dev scripts")
def test_run_sumcorr_dev_placeholder():
    pass
