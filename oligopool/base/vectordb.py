import sys

import time   as tt
import shutil as sh

from functools import wraps

from ShareDB import ShareDB

from . import utils as ut


class vectorDB:
    '''
    LevelDB based scalable on disk vectorDB for k-mer storage.
    '''

    def __init__(self, path, maximum_repeat_length):
        '''
        vectorDB constructor.

        :: path
           type - string
           desc - path to store vectorDB instance
        :: maximum_repeat_length
           type - integer
           desc - maximum shared repeat length

        Note:
            If reopening a vectorDB instance for reuse
            the maximum_repeat_length parameter is
            ignored, and the maximum_repeat_length
            from the vectorDB instance is used.
        '''

        # Aliasing
        maxreplen = maximum_repeat_length

        # Path Setup
        self.PATH = ut.removestarfix(
            string=path,
            fix='/',
            loc=1) + '/'
        if not self.PATH.endswith('vectorDB.ShareDB'):
            self.PATH += 'vectorDB.ShareDB'

        # Create/Open LevelDB object
        self.DB = ShareDB(
            path=self.PATH,
            map_size=None)

        # Length setup
        try:
            self.LEN = self.DB['LEN']
        except:
            self.LEN = 0
            self.DB['LEN'] = 0

        # K setup
        try:
            self.K = self.DB['K']
        except:
            self.K = int(maxreplen+1)
            self.DB['K'] = self.K

        # Verbosity Setup
        self.VERB = False

        # Object ALIVE status
        self.ALIVE = True

    def __repr__(self):
        '''
        User function to return the string representation of
        vectorDB.
        '''
        return 'vectorDB stored at {} with {} {}-mers'.format(
            self.PATH.removesuffix('vectorDB.ShareDB'), self.LEN, self.K)

    def __str__(self):
        '''
        See __repr__
        '''
        return self.__repr__()

    def alivemethod(method):
        '''
        Internal decorator to gate vectorDB operation
        once dropped. Internal use only.

        :: method
           type - function
           desc - function to gate
        '''
        @wraps(method)
        def wrapper(self, *args, **kwargs):
            if self.ALIVE:
                return method(self, *args, **kwargs)
            else:
                raise RuntimeError(
                    'vectorDB was closed or dropped')
        return wrapper

    def _verb_action(self, action, index, seq):
        '''
        Internal helper function to print sequence
        updates to vectorDB based on verbosity.
        Internal use only.

        :: action
           type - string
           desc - operation name
        :: index
           type - integer
           desc - sequence index name
        '''
        if self.VERB:
            psl = min(len(seq), 10)
            printed = ' {} Seq {}: {}...'.format(
                action, index, seq[:psl])
            clear_length = len(printed)
            sys.stdout.write(' '*clear_length+'\r')
            sys.stdout.write(printed)
            sys.stdout.flush()

    def _get(self, key):
        '''
        Internal helper to fetch value for given
        key, from vectorDB instance.
        Internal use only.

        :: key
           type - string
           desc - key to fetch from instance
        '''
        val = self.DB.get(key)
        if val is None:
            return None
        return val

    @alivemethod
    def add(self, seq, rna=True):
        '''
        User function to add seq k-mers to vectorDB.

        :: seq
           type - string
           desc - sequence to add to instance
        :: rna
           type - boolean
           desc - if True will convert seq to
                  DNA string prior to storage
        '''
        try:
            if not seq in ['K', 'LEN']:
                if rna:
                    seq = seq.replace('U', 'T')
                for kmer in ut.stream_canon_spectrum(
                    seq=seq, k=self.K):
                    if self._get(kmer) is None:
                        self.DB[kmer] = True
                        self.LEN += 1
            self.DB['LEN'] = self.LEN
        except Exception as E:
            raise E

    @alivemethod
    def multiadd(self, seq_list, rna=True):
        '''
        User function to add seq k-mers for each seq in
        seq_list to vectorDB.

        :: seq_list
           type - list
           desc - list of sequences to add to instance
        :: rna
           type - boolean
           desc - if True will convert seq to
                  DNA string prior to storage
        '''
        try:
            t0 = tt.time()
            pt = False
            index = 0
            for seq in seq_list:
                if not pt and self.VERB:
                    print('\n[Background Processing]')
                    pt = True
                if not seq in ['K', 'LEN']:
                    self._verb_action(
                        action='Adding',
                        index=index,
                        seq=seq)
                    index += 1
                    if rna:
                        seq = seq.replace('U', 'T')
                    for kmer in ut.stream_canon_spectrum(
                        seq=seq, k=self.K):
                        if self._get(kmer) is None:
                            self.DB[kmer] = True
                            self.LEN += 1
            self.DB['LEN'] = self.LEN
            self.DB.sync()
            if self.VERB:
                print('\n Time Elapsed: {:.2f} sec'.format(tt.time()-t0))
        except Exception as E:
            raise E

    @alivemethod
    def __contains__(self, seq, rna=True):
        '''
        Python dunder function to check existence of
        any k-mer from seq in vectorDB.

        :: seq
           type - string
           desc - sequence to check for existence
                  inside instance
        :: rna
           type - boolean
           desc - if True will convert seq to
                  DNA string prior to checking
        '''
        try:
            if rna:
                seq = seq.replace('U', 'T')
            for kmer in ut.stream_canon_spectrum(
                seq=seq, k=self.K):
                if self._get(kmer):
                    return True
            return False
        except Exception as E:
            raise E

    @alivemethod
    def multicheck(self, seq_list, rna=True):
        '''
        User function to check existence of any k-mer for
        each seq in seq_list in vectorDB.

        :: seq_list
           type - list
           desc - list of sequences to check for existence
                  inside instance
        :: rna
           type - boolean
           desc - if True will convert seq to
                  DNA string prior to checking
        '''
        try:
            for seq in seq_list:
                if rna:
                    seq = seq.replace('U', 'T')
                yield self.__contains__(
                    seq=seq,
                    rna=False)
        except Exception as E:
            raise E

    @alivemethod
    def __iter__(self):
        '''
        User fuction to iterate over k-mers stored in
        vectorDB.
        '''
        for key, _ in self.DB.items():
            if not key in ['K', 'LEN']:
                yield key

    @alivemethod
    def __len__(self):
        '''
        User function to return the number of keys stored
        in vectorDB.
        '''
        return self.LEN

    @alivemethod
    def remove(self, seq, rna=True):
        '''
        User function to remove all k-mers in seq
        from vectorDB.

        :: seq
           type - string
           desc - sequence to remove from instance
        :: rna
           type - boolean
           desc - if True will convert seq to
                  DNA string prior to removal
        '''
        try:
            if not seq in ['K', 'LEN']:
                if rna:
                    seq = seq.replace('U', 'T')
                for kmer in ut.stream_canon_spectrum(
                    seq=seq, k=self.K):
                    if self._get(kmer):
                        self.DB.remove(kmer)
                        self.LEN -= 1
                self.DB['LEN'] = self.LEN
        except Exception as E:
            raise E

    @alivemethod
    def multiremove(self, seq_list, rna=True, clear=False):
        '''
        User function to remove all k-mers from each
        seq in seq_list from vectorDB.

        :: seq
           type - list
           desc - list of sequences to remove
                  from instance
        :: rna
           type - boolean
           desc - if True will convert seq to
                  DNA string prior to removal
        '''
        try:
            t0 = tt.time()
            pt = False
            index = 0
            for seq in seq_list:
                if not pt and self.VERB:
                    print('\n[Background Processing]')
                    pt = True
                if not seq in ['K', 'LEN']:
                    self._verb_action(
                        action='Removing',
                        index=index,
                        seq=seq)
                    index += 1
                    if rna:
                        seq = seq.replace('U', 'T')
                    for kmer in ut.stream_canon_spectrum(
                        seq=seq, k=self.K):
                        if clear:
                            self.DB.remove(kmer)
                            self.LEN -= 1
                        else:
                            if self._get(kmer):
                                self.DB.remove(kmer)
                                self.LEN -= 1
            self.DB.sync()
            self.DB['LEN'] = self.LEN
            if self.VERB:
                print('\n Time Elapsed: {:.2f} sec'.format(tt.time()-t0))
        except Exception as E:
            raise E

    def clear(self, maxreplen=None):
        '''
        User function to clear all k-mers stored in
        vectorDB.

        :: maxreplen
           type - integer
           desc - a new maximum shared repeat length
                  to consider for future operations
        '''
        self.drop()
        self.DB = ShareDB(path=self.PATH, map_size=None)
        self.LEN = 0
        self.DB['LEN'] = 0
        if maxreplen is None:
            self.DB['K'] = self.K
        else:
            if not maxreplen is None:
                if not isinstance(maxreplen, int):
                    print('\n [ERROR]    maxreplen must be an integer, not {}'.format(
                        type(maxreplen)))
                    print(' [SOLUTION] Try correcting maxreplen\n')
                    raise ValueError
                self.K = int(maxreplen) + 1
                if maxreplen < 5:
                    print('\n [ERROR]    maxreplen must be greater than 4, not {}'.format(
                        maxreplen))
                    print(' [SOLUTION] Try correcting maxreplen\n')
                    raise ValueError
            self.DB['K'] = self.K
        self.ALIVE = True

    def close(self):
        '''
        User function to close vectorDB.
        '''
        if self.ALIVE:
            del self.DB
            self.DB = None
            self.ALIVE = False
            return True
        return False

    def drop(self):
        '''
        User function to drop vectorDB.
        '''
        if self.ALIVE:
            del self.DB
            self.DB = None
            sh.rmtree(self.PATH.removesuffix('vectorDB.ShareDB'))
            self.ALIVE = False
            return True
        return False