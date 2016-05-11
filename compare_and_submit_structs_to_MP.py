import pymongo
from pymatgen import MPRester, Composition, Structure
from pymatgen.matproj.snl import StructureNL

client = pymongo.MongoClient()
db = client.springer
mpr = MPRester()


def create_mincoll():
    origcoll = db['pauling_file']
    min_collname = 'pauling_file_mpmin'
    db[min_collname].drop()
    origcoll.aggregate([{'$match': {'structure': {'$exists': True}, 'metadata._structure.is_ordered': True,
                                    'metadata._structure.is_valid': True, 'errors': {
            '$in': ['structural composition and refined/alphabetic formula do not match']}}},
                        {'$project': {'key': 1, 'metadata': 1, 'structure': 1, 'webpage_link': 1}},
                        {'$out': min_collname}])
    db[min_collname].create_index([('key', pymongo.ASCENDING)], unique=True)
    # Remove Deuterium
    for doc in db[min_collname].find().batch_size(75):
        for el in doc['metadata']['_structure']['elements']:
            if el == 'D':
                db[min_collname].remove({'key': doc['key']})
                break


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
    create_mincoll()
    mp_unique_comps = set()
    pf_unique_comps = set()
    coll = db['pauling_file_mpmin']
    new_coll = db['pf_to_mp']
    new_coll.drop()
    mp_comps = mpr.query(criteria={}, properties=["snl_final.reduced_cell_formula_abc"])
    print 'Total number of MP comps = {}'.format(len(mp_comps))
    for mp_comp in mp_comps:
        mp_unique_comps.add(mp_comp["snl_final.reduced_cell_formula_abc"])
    print 'Number of unique MP comps = {}'.format(len(mp_unique_comps))
    x = 0
    for doc in coll.find().batch_size(75):
        x += 1
        if x % 1000 == 0:
            print x
        # if x > 1:
        #     break
        pf_unique_comps.add(doc['metadata']['_structure']['reduced_cell_formula_abc'])
        if doc['metadata']['_structure']['reduced_cell_formula_abc'] not in mp_unique_comps:
            new_coll.insert(doc_to_snl(doc).as_dict())
            # with open('PaulingFile_example.json', 'w') as outfile:
            #     json.dump(doc_to_snl(doc).as_dict(), outfile)
    new_coll.create_index([('about._pauling_file.key', pymongo.ASCENDING)], unique=True)
    print 'Number of PF unique comps = {}'.format(len(pf_unique_comps))
    new_comps = pf_unique_comps.difference(mp_unique_comps)
    print 'Number of new compositions in PF = {}'.format(len(new_comps))
    '''
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
