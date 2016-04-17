import pymongo
import json
from pymatgen import MPRester, Composition, Structure
from pymatgen.matproj.snl import StructureNL

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


def doc_to_snl(doc):
    return StructureNL(Structure.from_dict(doc['structure']), [{"name": "Saurabh Bajaj", "email": "sbajaj@lbl.gov"},
                                                               {"name": "Anubhav Jain", "email": "ajain@lbl.gov"}], [],
                       '', ['Pauling file'], {'_pauling_file': {'key': doc['key'],
                                                                'reference': doc['metadata']['_Springer']['geninfo'][
                                                                    'ref'], 'general_info': {
                'Sample Detail(s)': doc['metadata']['_Springer']['geninfo']['Sample Detail(s)'],
                'Standard Formula': doc['metadata']['_Springer']['geninfo']['Standard Formula'],
                'Mineral Name(s)': doc['metadata']['_Springer']['geninfo']['Mineral Name(s)'],
                'Phase Prototype': doc['metadata']['_Springer']['geninfo']['Phase Prototype'],
                'Structure Class(es)': doc['metadata']['_Springer']['geninfo']['Structure Class(es)'],
                'Measurement Detail(s)': doc['metadata']['_Springer']['geninfo']['Measurement Detail(s)'],
                'Phase Label(s)': doc['metadata']['_Springer']['geninfo']['Phase Label(s)']},
                                                                'expdetails': doc['metadata']['_Springer'][
                                                                    'expdetails'],
                                                                'title': doc['metadata']['_Springer']['title']}},
                       [{'name': 'Pauling file', 'url': doc['webpage_link'], 'description': {'key': doc['key']}}])


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
    pf_unique_comps = set()
    pf_unique_comps_structs = []
    mp_comps = set()
    mp_unique_comps = set()
    coll = db['pauling_file']
    x = 0
    for doc in coll.find({'structure': {'$exists': True}}).batch_size(75):
        if doc['metadata']['_structure']['is_ordered']:
            x += 1
            # if x % 1000 == 0:
            #     print x
            if x > 1:
                break
            pf_unique_comps.add(doc['metadata']['_structure']['reduced_cell_formula_abc'])
            print doc_to_snl(doc).as_dict()
            with open('PaulingFile_example.json', 'w') as outfile:
                json.dump(doc_to_snl(doc).as_dict(), outfile)
    print 'Number of PF unique comps = {}'.format(len(pf_unique_comps))
    '''
    mp_comps = mpr.query(criteria={}, properties=["pretty_formula"])
    print 'Number of MP comps = {}'.format(len(mp_comps))
    for mp_comp in mp_comps:
        mp_unique_comps.add(Composition(mp_comp['pretty_formula']).alphabetical_formula)
    print 'Number of MP unique comps = {}'.format(len(mp_unique_comps))
    new_comps = pf_unique_comps.difference(mp_unique_comps)
    print 'Number of new compositions in PF = {}'.format(len(new_comps))
    for s in structures:
        found = mpr.find_structure(s)
        print found
        if len(found) > 0:
            mp_structs.extend(found)
        else:
            new_structures.append(s)
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
    '''
