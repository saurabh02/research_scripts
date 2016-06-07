import os
from pymatgen import MPRester
from collections import defaultdict
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
from sklearn.metrics import mean_squared_error
import pickle
import numpy as np

mpr = MPRester()

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(module_dir, 'data_bondlengths')


def get_bondlengths(structure_lst):
    bond_lengths = defaultdict(list)
    for struct in structure_lst:
        for site_idx, site in enumerate(struct.sites):
            try:
                voronoi_sites = VoronoiCoordFinder(struct).get_coordinated_sites(site_idx)
            except RuntimeError as r:
                print 'Error for site {} in {}: {}'.format(site, struct.composition, r)
                continue
            except ValueError as v:
                print 'Error for site {} in {}: {}'.format(site, struct.composition, v)
                continue
            for vsite in voronoi_sites:
                # site_dist = site.distance(vsite)
                # print site_dist
                s_dist = np.linalg.norm(vsite.coords - site.coords)
                if s_dist < 0.1:
                    continue
                bond = '-'.join(sorted([site.species_string, vsite.species_string]))
                bond_lengths[bond].append(s_dist)
    return bond_lengths


def get_avg_bondlengths(bondlengths):
    avg_bondlengths = {}
    for bond in bondlengths:
        avg_bondlengths[bond] = sum(bondlengths[bond])/len(bondlengths[bond])
    return avg_bondlengths


def save_avg_bondlengths(nelements, e_above_hull=0.05):
    criteria = {'nelements': {'$lte': nelements}}
    mp_results = mpr.query(criteria=criteria, properties=['task_id', 'e_above_hull', 'structure'])
    mp_structs = []
    for i in mp_results:
        if i['e_above_hull'] < e_above_hull:
            mp_structs.append(i['structure'])
    bls = get_bondlengths(mp_structs)
    avg_bls = get_avg_bondlengths(bls)
    with open(os.path.join(data_dir, 'nelements_' + str(nelements) + '_avgbls.pkl'), 'w') as f:
        pickle.dump(avg_bls, f, pickle.HIGHEST_PROTOCOL)


def get_rmse(avg_bls_db, bond_lengths):
    rmse = 0
    for bond in bond_lengths:
        rmse += mean_squared_error(bond_lengths[bond], [avg_bls_db[bond]]*len(bond_lengths[bond]))**0.5
    return rmse


def predict_volume(structure):
    with open(os.path.join(data_dir, 'nelements_2_avgbls.pkl'), 'r') as f:
        avg_bls_db = pickle.load(f)
    starting_volume = structure.volume
    predicted_volume = starting_volume
    min_rmse = get_rmse(avg_bls_db, get_bondlengths([structure]))
    for i in range(81, 121):
        test_volume = (i * 0.01) * starting_volume
        structure.scale_lattice(test_volume)
        test_structure_bls = get_bondlengths([structure])
        test_rmse = get_rmse(avg_bls_db, test_structure_bls)
        if test_rmse < min_rmse:
            min_rmse = test_rmse
            predicted_volume = test_volume
    return predicted_volume, min_rmse


if __name__ == '__main__':
    # save_avg_bondlengths(2)
    new_fe_struct = mpr.get_structure_by_material_id('mp-30434')
    print 'Starting volume = {}'.format(new_fe_struct.volume)
    print 'Predicted volume = {} with RMSE = {}'.format(predict_volume(new_fe_struct)[0], predict_volume(new_fe_struct)[1])

