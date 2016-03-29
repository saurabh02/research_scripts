import pymongo
from pymatgen import Structure, MPRester, Composition
from pymatgen.matproj.snl import StructureNL
import pickle

client = pymongo.MongoClient()
db = client.springer
mpr = MPRester()


def get_meta_from_structure(structure):
    """
    Used by `structure_to_mock_job`, to "fill out" a job document.
    :param structure: pymatgen structure object
    :return: (dict) structure metadata
    """
    comp = structure.composition
    elsyms = sorted(set([e.symbol for e in comp.elements]))
    meta = {'nsites': len(structure),
            'elements': elsyms,
            'nelements': len(elsyms),
            'formula': comp.formula,
            'reduced_cell_formula': comp.reduced_formula,
            'reduced_cell_formula_abc': Composition(comp.reduced_formula).alphabetical_formula,
            'anonymized_formula': comp.anonymized_formula,
            'chemsystem': '-'.join(elsyms),
            'is_ordered': structure.is_ordered,
            'is_valid': structure.is_valid()}
    return meta


def structure_to_mock_job(structure):
    # Needs at least one author. This is for a mock job, so can put whatever.
    snl = StructureNL(structure, [{"name": "Saurabh Bajaj", "email": "sbajaj@lbl.gov"},
                                  {"name": "Anubhav Jain", "email": "ajain@lbl.gov"}])
    job = snl.as_dict()
    if 'is_valid' not in job:
        job.update(get_meta_from_structure(snl.structure))
    sorted_structure = snl.structure.get_sorted_structure()
    job.update(sorted_structure.as_dict())
    return job


def job_is_submittable(job):
    snl = StructureNL.from_dict(job)
    # mpworks.processors.process_submissions.SubmissionProcessor#submit_new_workflow
    max_sites = 200  # SubmissionProcessor.MAX_SITES above
    # from mpworks.workflows.wf_utils import NO_POTCARS
    no_potcars = ['Po', 'At', 'Rn', 'Fr', 'Ra', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
    if len(snl.structure.sites) > max_sites:
        print 'REJECTED WORKFLOW FOR {} - too many sites ({})'.format(
            snl.structure.formula, len(snl.structure.sites))
    elif not job['is_valid']:
        print 'REJECTED WORKFLOW FOR {} - invalid structure (atoms too close)'.format(
            snl.structure.formula)
    elif len(set(no_potcars) & set(job['elements'])) > 0:
        print 'REJECTED WORKFLOW FOR {} - invalid element (No POTCAR)'.format(
            snl.structure.formula)
    elif not job['is_ordered']:
        print 'REJECTED WORKFLOW FOR {} - invalid structure (disordered)'.format(
            snl.structure.formula)
    else:
        return True
    return False


if __name__ == '__main__':
    structures = []
    mp_structs = []
    new_structures = []
    coll = db['pauling_file']
    x = 0
    for doc in coll.find({'structure': {'$exists': True}}).batch_size(75):
        x += 1
        if x % 1000 == 0:
            print x
        structures.extend([Structure.from_dict(doc['structure'])])
    print 'Number of extracted structures = {}'.format(len(structures))
    with open('PF_structures', 'w') as w:
        pickle.dump(structures, w)
    num_found = 0
    num_new = 0
    for s in structures:
        found = mpr.find_structure(s)
        if len(found) > 0:
            mp_structs.extend(found)
            num_found += 1
            if num_found % 100 == 0:
                print num_found
        else:
            new_structures.append(s)
            num_new += 1
            if num_new % 100 == 0:
                print num_new
    if len(mp_structs) > 0:
        print("Number of filtered out structures already on MP: {}".format(len(mp_structs)))
    if len(new_structures) > 0:
        print("Number of new structures: {}".format(len(mp_structs)))
    with open('mp_structs', 'w') as w:
        pickle.dump(mp_structs, w)
    with open('new_structs', 'w') as w:
        pickle.dump(new_structures, w)
    num_subm = 0
    submittables = []
    for s in new_structures:
        if job_is_submittable(structure_to_mock_job(s)):
            submittables.append(s)
            num_subm += 1
            if num_subm % 100 == 0:
                print num_subm
    print("Number of new submittable structures: {}".format(len(submittables)))
    with open('submittables', 'w') as w:
        pickle.dump(submittables, w)

