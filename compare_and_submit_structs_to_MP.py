import pymongo
from pymatgen import Structure, MPRester

client = pymongo.MongoClient()
db = client.springer
mpr = MPRester()


if __name__ == '__main__':
    structures = []
    mp_ids = []
    new_structures = []
    coll = db['pauling_file']
    x = 0
    for doc in coll.find({'structure': {'$exists': True}}).batch_size(75):
        x += 1
        if x % 1000 == 0:
            print x
        structures.extend(Structure.from_dict(doc['structure']))
    print 'Number of extracted structures = {}'.format(len(structures))
    for s in structures:
        found = mpr.find_structure(s)
        if len(found) > 0:
            mp_ids.extend(found)
        else:
            new_structures.append(s)
    if len(mp_ids) > 0:
        print("Number of filtered out structures already on MP: {}".format(len(mp_ids)))
    if len(new_structures) > 0:
        print("Number of new structures: {}".format(len(mp_ids)))