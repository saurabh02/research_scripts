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


class VolumePredictor(object):
    """
    Class to predict the volume of a given structure, based on averages of bond lengths of given input structures.
    """

    def __init__(self):
        self.bond_lengths = defaultdict(list)
        self.avg_bondlengths = {}

    def get_bondlengths(self, structure):
        """
        Get all bond lengths in a structure by bond type.

        :param structure: pymatgen structure object
        :return: (defaultdict(list)) with bond types in the format 'A-B' as keys, and bond lengths as values
        """
        bondlengths = defaultdict(list)
        for site_idx, site in enumerate(structure.sites):
            try:
                voronoi_sites = VoronoiCoordFinder(structure).get_coordinated_sites(site_idx)
            except RuntimeError as r:
                print 'Error for site {} in {}: {}'.format(site, structure.composition, r)
                continue
            except ValueError as v:
                print 'Error for site {} in {}: {}'.format(site, structure.composition, v)
                continue
            for vsite in voronoi_sites:
                s_dist = np.linalg.norm(vsite.coords - site.coords)
                if s_dist < 0.1:
                    continue
                bond = '-'.join(sorted([site.species_string, vsite.species_string]))
                bondlengths[bond].append(s_dist)
        return bondlengths

    def fit(self, structures, volumes):
        """
        Given a set of input structures, it stores bond lengths in a defaultdict(list) and average bond lengths in a
        dictionary by bond type (keys).

        :param structures: (list) list of pymatgen structure objects
        :param volumes: (list) corresponding list of volumes of structures
        """
        for struct in structures:
            struct_bls = self.get_bondlengths(struct)
            for bond in struct_bls:
                self.bond_lengths[bond].extend(struct_bls[bond])
        for bond in self.bond_lengths:
            # print self.bond_lengths[bond]
            self.avg_bondlengths[bond] = sum(self.bond_lengths[bond])/len(self.bond_lengths[bond])

    def get_rmse(self, structure_bls):
        """
        Calculates root mean square error between bond lengths in a structure and the average expected bond lengths.

        :param structure_bls: defaultdict(list) of bond lengths as output by the function "get_bondlengths()"
        :return: (float) root mean square error
        """
        rmse = 0
        for bond in structure_bls:
            rmse += mean_squared_error(structure_bls[bond], [self.avg_bondlengths[bond]]*len(structure_bls[bond]))**0.5
        return rmse

    def predict(self, structure):
        """
        Predict volume of a given structure based on rmse against average bond lengths from a set of input structures.
        Note: run this function after initializing the "self.avg_bondlengths" variable, through either of the functions
        fit() or get_avg_bondlengths().

        :param structure: (structure) pymatgen structure object to predict the volume of
        :return: (tuple) predited volume (float) and its rmse (float)
        """
        starting_volume = structure.volume
        predicted_volume = starting_volume
        min_rmse = self.get_rmse(self.get_bondlengths(structure))
        for i in range(81, 121):
            test_volume = (i * 0.01) * starting_volume
            structure.scale_lattice(test_volume)
            test_rmse = self.get_rmse(self.get_bondlengths(structure))
            if test_rmse < min_rmse:
                min_rmse = test_rmse
                predicted_volume = test_volume
        return predicted_volume, min_rmse

    def save_avg_bondlengths(self, filename):
        """
        Save the average bond lengths calculated by fit().

        :param filename: name of file to store to
        """
        with open(os.path.join(data_dir, filename), 'w') as f:
            pickle.dump(self.avg_bondlengths, f, pickle.HIGHEST_PROTOCOL)

    def get_avg_bondlengths(self, filename):
        """
        Extract the saved average bond lengths and save them in the class variable "self.avg_bondlengths".

        :param filename: name of file to extract average bond lengths from
        :return:
        """
        with open(os.path.join(data_dir, filename), 'r') as f:
            self.avg_bondlengths = pickle.load(f)


def get_bondlengths(structure_lst):
    bondlengths = defaultdict(list)
    for structure in structure_lst:
        for site_idx, site in enumerate(structure.sites):
            try:
                voronoi_sites = VoronoiCoordFinder(structure).get_coordinated_sites(site_idx)
            except RuntimeError as r:
                print 'Error for site {} in {}: {}'.format(site, structure.composition, r)
                continue
            except ValueError as v:
                print 'Error for site {} in {}: {}'.format(site, structure.composition, v)
                continue
            for vsite in voronoi_sites:
                s_dist = np.linalg.norm(vsite.coords - site.coords)
                if s_dist < 0.1:
                    continue
                bond = '-'.join(sorted([site.species_string, vsite.species_string]))
                bondlengths[bond].append(s_dist)
    return bondlengths


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
    new_struct = mpr.get_structure_by_material_id('mp-97')
    print 'Starting volume = {}'.format(new_struct.volume)
    pred_vol = predict_volume(new_struct)
    print 'Predicted volume = {} with RMSE = {}'.format(pred_vol[0], pred_vol[1])
    criteria = {'nelements': {'$lte': 2}}
    mp_results = mpr.query(criteria=criteria, properties=['task_id', 'e_above_hull', 'structure'])
    mp_structs = []
    mp_vols = []
    for i in mp_results:
        if i['e_above_hull'] < 0.05:
            mp_structs.append(i['structure'])
            mp_vols.append(i['structure'].volume)
    pv = VolumePredictor()
    pv.fit(mp_structs, mp_vols)
    pv.save_avg_bondlengths()
    # a = pv.predict(new_struct)
    # print 'Predicted volume = {} with RMSE = {}'.format(a[0], a[1])


