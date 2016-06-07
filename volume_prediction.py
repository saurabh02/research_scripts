import os
from pymatgen import MPRester
from collections import defaultdict
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
from sklearn.metrics import mean_squared_error
import pickle

mpr = MPRester()

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(module_dir, 'data_bondlengths')


def get_bondlengths(structure_lst, radius=4):
    bond_lengths = defaultdict(list)
    for struct in structure_lst:
        for site_idx, site in enumerate(struct.sites):
            voronoi_sites = VoronoiCoordFinder(struct, cutoff=4).get_coordinated_sites(site_idx)
            for vsite in voronoi_sites:
                site_dist = site.distance(vsite)
                if site_dist < 0.1:
                    continue
                bond = '-'.join(sorted([site.species_string, vsite.species_string]))
                bond_lengths[bond].append(site_dist)
    return bond_lengths


def get_avg_bondlengths(bondlengths):
    avg_bondlengths = {}
    for bond in bondlengths:
        avg_bondlengths[bond] = sum(bondlengths[bond])/len(bondlengths[bond])
    return avg_bondlengths


def save_avg_bondlengths(nelements, e_above_hull=0.05):
    criteria = {'nelements': {'$lte': nelements}}
    mp_results = mpr.query(criteria=criteria, properties=['task_id', 'e_above_hull'])
    mp_results_stable = []
    for i in mp_results:
        if i['e_above_hull'] < e_above_hull:
            struc_result = mpr.query({'task_id': i['task_id']}, ['structure'])
            mp_results_stable.append(struc_result[0]['structure'])
    bls = get_bondlengths(mp_results_stable)
    avg_bls = get_avg_bondlengths(bls)
    with open(os.path.join(data_dir, 'nelements_' + str(nelements) + '_avgbls.pkl'), 'wb') as f:
        pickle.dump(avg_bls, f, pickle.HIGHEST_PROTOCOL)


def get_rmse(avg_bls_db, bond_lengths):
    rmse = {}
    for bond in bond_lengths:
        rmse[bond] = mean_squared_error(bond_lengths[bond], [avg_bls_db[bond]]*len(bond_lengths[bond]))**0.5
    return rmse


def predict_volume(structure):
    with open(os.path.join(data_dir, 'nelements_2_avgbls.pkl'), 'r') as f:
        avg_bls_db = pickle.load(f)
    starting_volume = structure.volume
    predicted_volume = starting_volume
    curr_min_rmse = get_rmse(avg_bls_db, get_bondlengths([structure])).values()[0]      # TODO: need to account for more than 1 bond type
    print curr_min_rmse
    for i in range(81, 121):
        test_volume = (i * 0.01) * starting_volume
        structure.scale_lattice(test_volume)
        test_structure_bls = get_bondlengths([structure])
        test_rmse = get_rmse(test_structure_bls).values()[0]       # TODO: need to account for more than 1 bond type
        if test_rmse < curr_min_rmse:
            curr_min_rmse = test_rmse
            predicted_volume = test_volume
    return predicted_volume


if __name__ == '__main__':
    # save_avg_bondlengths(2)
    # with open(os.path.join(data_dir, 'nelements_' + str(2) + '_avgbls.pkl'), 'rb') as f:
    #     avg_bls = pickle.load(f)
    # print avg_bls.keys()
    # '''
    new_fe_struct = mpr.query(criteria={'task_id': 'mp-568345'},
                           properties=['task_id', 'pretty_formula', 'structure', 'e_above_hull'])[0]['structure']
    print get_bondlengths([new_fe_struct])
    # print 'Starting volume = {}'.format(new_fe_struct.volume)
    # print 'Predicted volume = {}'.format(predict_volume(new_fe_struct))
    # '''

