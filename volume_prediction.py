from pymatgen import MPRester
from collections import defaultdict
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
from sklearn.metrics import mean_squared_error
import pickle

mpr = MPRester()


def get_bondlengths(structure_lst):
    bond_lengths = defaultdict(list)
    for struct in structure_lst:
        for site_idx, site in enumerate(struct.sites):
            # voronoi_polyhedra = VoronoiCoordFinder(structure, cutoff=4).get_voronoi_polyhedra(site_idx)
            # print 'Voronoi polyhedra for site {} = {}'.format(site, voronoi_polyhedra)
            # for ngbsite in voronoi_polyhedra:
                # print structure.get_distance(site_idx, ngbsite)
                # print ngbsite
            sites_in_sphere = struct.get_sites_in_sphere(site.coords, r=4, include_index=True)
            for s in sites_in_sphere:
                if s[1] == 0:
                    continue
                bond = '-'.join(sorted([site.species_string, s[0].species_string]))
                # print '{} = {}'.format(bond, s[1])
                bond_lengths[bond].append(s[1])
        # print bond_lengths
        # print get_rmse(bond_lengths)
    return bond_lengths


def get_avg_bondlengths(bondlengths):
    avg_bondlengths = {}
    for bond in bondlengths:
        avg_bondlengths[bond] = sum(bondlengths[bond])/len(bondlengths[bond])
    return avg_bondlengths


def get_rmse(bond_lengths):
    with open('element_avg_bl.pkl', 'rb') as f:
        avg_bls = pickle.load(f)
    rmse = {}
    for bond in bond_lengths:
        rmse[bond] = mean_squared_error(bond_lengths[bond], [avg_bls[bond]]*len(bond_lengths[bond]))**0.5
    return rmse


def predict_volume(structure):
    starting_volume = structure.volume
    predicted_volume = starting_volume
    starting_bls = get_bondlengths([structure])
    curr_min_rmse = get_rmse(starting_bls).values()[0]      # TODO: need to account for more than 1 bond type
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
    # criteria = {'nelements': {'$in': [1,2]}}
    # criteria = {'task_id': {'$in': ['mp-656887', 'mp-1703']}}
    # criteria = {'nelements': {'$in': [1]}, 'elements': {'$in': ['Fe']}}
    '''
    criteria = {'nelements': {'$in': [1]}}
    mp_results = mpr.query(criteria=criteria,
                           properties=['task_id', 'pretty_formula', 'structure', 'e_above_hull'])
    mp_results_stable = []
    mp_results_stablestructs = []
    for i in mp_results:
        if i['e_above_hull'] < 0.5:
            mp_results_stable.append(i)
            mp_results_stablestructs.append(i['structure'])
    bl = get_bondlengths(mp_results_stablestructs)
    avg_bl = get_avg_bondlengths(bl)
    with open('element_avg_bl.pkl', 'wb') as f:
        pickle.dump(avg_bl, f, pickle.HIGHEST_PROTOCOL)
    '''
    new_fe_struct = mpr.query(criteria={'task_id': 'mp-568345'},
                           properties=['task_id', 'pretty_formula', 'structure', 'e_above_hull'])[0]['structure']
    print 'Starting volume = {}'.format(new_fe_struct.volume)
    print 'Predicted volume = {}'.format(predict_volume(new_fe_struct))
    # for i in mp_results_stable:
        # get_bond_lenghts_in_structure(i['structure'])
        # bl = get_bondlengths(i['structure'])

