import os

from . import utils as ut
from . import vectordb as db


def normalize_background_dir(path):
    '''
    Normalize a background directory path,
    stripping vectorDB.ShareDB suffix.
    Internal use only.

    :: path
       type - string
       desc - background directory path
    '''
    p = path.rstrip('/\\')
    if p.endswith('vectorDB.ShareDB'):
        p = os.path.dirname(p)
    if 'vectorDB.ShareDB' in p:
        head, _sep, _tail = p.partition('vectorDB.ShareDB')
        p = head.rstrip('/\\')
    return p

def _sharedb_dir(path):
    '''
    Compute the ShareDB directory for a background DB.
    Internal use only.

    :: path
       type - string
       desc - background directory path
    '''
    root = normalize_background_dir(path)
    return os.path.join(root, 'vectorDB.ShareDB')

def is_lmdb_dir(path):
    '''
    Check if path is an LMDB environment directory.
    Internal use only.

    :: path
       type - string
       desc - candidate LMDB directory
    '''
    if not os.path.isdir(path):
        return False
    # Fast path: file presence check.
    if not os.path.isfile(os.path.join(path, 'data.mdb')):
        return False
    if not os.path.isfile(os.path.join(path, 'ShareDB.config')):
        return False
    # Strong check: attempt a readonly open (no creation).
    try:
        import lmdb
        env = lmdb.open(
            path,
            subdir=True,
            readonly=True,
            create=False,
            lock=False,
            readahead=False,
            max_readers=1,
            max_dbs=0)
        env.close()
        return True
    except Exception:
        return False

def is_vectordb_dir(path):
    '''
    Check if path contains a vectorDB.ShareDB LMDB directory.
    Internal use only.

    :: path
       type - string
       desc - background directory path
    '''
    return is_lmdb_dir(_sharedb_dir(path))

def inspect_background(path):
    '''
    Inspect a background DB directory.
    Internal use only.

    :: path
       type - string
       desc - background directory path

    Returns:
        (ok, warns, meta)
    '''
    warns = {}
    root = normalize_background_dir(path)
    if not os.path.isdir(root):
        return False, {}, {'path': root, 'reason': 'directory_missing'}
    if not is_vectordb_dir(root):
        return False, {}, {'path': root, 'reason': 'not_a_vectordb_directory'}
    try:
        # At this point we have confirmed LMDB looks valid; now use ShareDB-backed
        # vectorDB to load the canonical K/LEN metadata.
        vdb = db.vectorDB(path=root, maximum_repeat_length=0)
        meta = {'directory': root, 'K': int(vdb.K), 'LEN': int(vdb.LEN)}
        vdb.close()
        return True, warns, meta
    except Exception as exc:
        return False, warns, {'directory': root, 'reason': str(exc)}

def inspect_zip(path, load_dict_names):
    '''
    Inspect a zip archive produced by analysis modules.
    Internal use only.

    :: path
       type - string
       desc - archive file path
    :: load_dict_names
       type - tuple / list
       desc - archived dict names to attempt to load

    Returns:
        (ok, warns, meta)
    '''
    warns = {}
    meta = {'file': path, 'entries': []}
    try:
        archive = ut.get_archive(path)
    except Exception as exc:
        return False, warns, {'file': path, 'reason': 'failed_to_open_archive: {}'.format(exc)}
    try:
        info = archive.infolist()
        meta['entries'] = [{'name': i.filename, 'bytes': int(i.file_size)} for i in info]
        names = {i.filename for i in info}
        meta['schema'] = {'expected': list(load_dict_names), 'present': sorted(names.intersection(load_dict_names))}
        missing = sorted(set(load_dict_names).difference(names))
        if missing:
            warns['missing_entries'] = missing

        for name in load_dict_names:
            if name not in names:
                continue
            try:
                meta[name] = ut.loaddict(archive=archive, dfile=name)
            except Exception as exc:
                warns.setdefault('failed_load', {})[name] = str(exc)

        return True, warns, meta
    finally:
        archive.close()
