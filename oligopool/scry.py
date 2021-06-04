import collections as cx
import numpy as np
import numba as nb
import edlib as ed


class Scry(object):

    '''
    "The Inner Eye does not See upon command!"
    - Professor Sybill Trelawney
    '''

    def __init__(self):
        '''
        Scry classifier constructor.
        Internal use only.
        '''

        # Maps
        self.c = None # Corpus Map
        self.X = None # Corpus List
        self.d = None # Spectrum Dictionary

        # Instance Variables
        self.n = None # Corpus Length
        self.k = None # Corpus k-value
        self.t = None # Corpus t-value
        self.z = None

        self.pt = None # Primed t-value
        self.pm = None # Primed mode

        # Sentinels
        self.trained = False
        self.primed  = False

    @staticmethod
    def absorb(map1, map2):
        '''
        Absorb map1 into map2.
        Internal use only.

        :: map1
           type - dict
           desc - source map
        :: map2
           type - dict
           desc - destination map
        '''

        while map1:
            k,v = map1.popitem()
            map2[k] = v
        return map2

    @staticmethod
    def stream_kmers(x, n, k, locations=False):
        '''
        Stream all k-mers from sequence x.
        Internal use only.

        :: x
           type - string
           desc - sequence string x
        :: n
           type - integer
           desc - length of x
        :: k
           type - integer
           desc - length of k-mers
        :: locations
           type - boolean
           desc - if True will stream k-mer
                  indices along with k-mers
                  (default=False)
        '''

        if locations:
            return ((x[i:i+k],i) for i in range(n-k+1))
        return (x[i:i+k] for i in range(n-k+1))

    @staticmethod
    def stream_reduced(store):
        '''
        Stream all substrings reduced by length
        of 1 from store. Internal use only.

        :: store
           type - list
           desc - contig storage
        '''

        yield store[0][:-1]
        yield store[0][+1:]
        i = 1
        while i < len(store):
            yield store[i][+1:]
            i += 1

    @staticmethod
    def stream_substrings(x, p, q):
        '''
        Stream all substrings of length
        p through q from sequence x.
        Internal use only.

        :: x
           type - string
           desc - sequence to stream
                  contigs from
        :: p
           type - integer
           desc - starting contig length
                  (inclusive)
        :: q
           type - integer
           desc - final contig length
                  (inclusive)
        '''

        # Initial Check
        if q > p:
            return None

        # Stream Initial Substrings
        store = []
        for kmer in Scry.stream_kmers(
            x=x, n=len(x), k=p):
            yield kmer, p
            store.append(kmer)

        # Stream Reduced Substrings
        l = p-1
        while l >= q:
            _store = []
            for reduced in Scry.stream_reduced(
                store=store):
                yield reduced, l
                _store.append(reduced)
            store = _store
            l -= 1

    @staticmethod
    def index(X, Y, n, k, t, liner=None):
        '''
        Index given sequences and labels
        in iterables. Internal use only.

        :: X
           type - iterable
           desc - iterable of all sequences
                  to index
        :: Y
           type - iterable
           desc - iterable of labels of all
                  indexed sequences
        :: k
           type - integer
           desc - length of k-mers
        :: t
           type - integer
           desc - allowed number of k-mer
                  location displacement
        :: liner
           type - coroutine / None
           desc - dynamic printing
        '''

        # Book-keeping
        c = {}
        d = {}
        B = set()

        # Indexing Phase 1 - Build Initial Index Objects
        for x,y in zip(X,Y):

            # Show Update
            if not liner is None:
                if y%100 == 0:
                    liner.send(
                        ' Barcode Model: Phase 1 on {}'.format(
                            'Barcode {:,}'.format(
                                y)))

            # Store Substrings in Corpus Map
            for ss,_ in Scry.stream_substrings(
                x=x, p=n, q=n-t):
                if not ss in c:
                    c[ss] = y
                else:
                    c[ss] = None # Duplicate Sentinel

            # Store k-mers in Spectrum Dictionary
            if t > 0: # Mutation Expected in Query
                for kmer,idx in Scry.stream_kmers(
                    x=x, n=n, k=k, locations=True):
                    if not kmer in d:
                        d[kmer] = []
                    d[kmer].append((idx,y))

        # Indexing Phase 2 - Optimize Spectrum Dictionary
        if t > 0: # Mutations Expected in Query
            for kcnt,kmer in enumerate(d):

                # Show Update
                if not liner is None:
                    liner.send(
                        ' Barcode Model: Phase 2 on {}'.format(
                            'k-mer {:,}'.format(
                                kcnt)))

                # Sort Current Record
                d[kmer].sort()

                # Process All Entries
                idxgroup = [None]*(n-k+1)
                ygroup   = []
                previdx  = None
                for ix,(idx,y) in enumerate(d[kmer]):
                    ygroup.append(y)
                    if idx != previdx:
                        idxgroup[idx] = [ix, ix+1]
                        previdx = idx
                    else:
                        idxgroup[idx][-1] += 1

                # Update k-mer Record
                d[kmer] = (np.array(ygroup),
                    tuple(tuple(r) if not r is None else r \
                          for r in idxgroup))

        # Return Objects
        return c,d

    def train(self, X, Y, n, k, t, liner=None):
        '''
        Train a Scry model.

        :: X
           type - iterable
           desc - iterable of all sequences
                  to fit
        :: Y
           type - iterable of labels of all
                  trained sequences
        :: n
           type - integer
           desc - length of each sequence
        :: k
           type - integer
           desc - k-value for sequences
        :: t
           type - integer
           desc - maximum mismatches in
                  sequence variants
        :: liner
           type - coroutine / None
           desc - dynamic printing
        '''

        X = list(X)
        Y = list(Y)
        c,d = self.index(
            X=X, Y=Y, n=n, k=k, t=t,
            liner=liner)
        self.X = X
        self.c = self.absorb(c, {})
        self.d = self.absorb(d, {})
        self.n = n
        self.k = k
        self.t = t
        self.z = len(X)
        self.trained = True
        return self

    def prime(self, t=0, mode=0):
        '''
        Prime instance for prediction.

        :: t
           type - integer
           desc - considered number of
                  mismatches in variants
                  (default=0)
        :: mode
           type - integer
           desc - prediction task mode
                  0 = fast / near-exact mode
                  1 = slow / sensitive  mode
                  (default=0)
        '''

        # Did we train this instance?
        if not self.trained:
            raise RuntimeError(
                'Untrained Scry instance')

        # Is mode valid?
        if (mode != 0) and (mode != 1):
            raise RuntimeError(
                'Invalid Scry prediction mode: {}'.format(
                    mode))

        # How many mismatches allowed?
        t = self.t if t > self.t else t if t >= 0 else 0

        # Store Primed Records
        self.pt = t
        self.pm = mode
        self.primed = True

    @staticmethod
    @nb.njit
    def update_scorevector(sv, iv, sr, en, sk):
        '''
        Update scorevector indices by k.
        Internal use only.

        :: sv
           type - np.array
           desc - score vector
        :: iv
           type - np.array
           desc - index vector
        :: sr
           type - integer
           desc - starting index (inclusive)
        :: en
           type - integer
           desc - ending index (inclusive)
        :: sk
           type - integer
           desc - score value
        '''
        for r in iv[sr:en]:
            sv[r] += sk

    @staticmethod
    @nb.njit
    def select_maxindex(sv):
        '''
        Select index with maximum score.
        Internal use only.

        :: sv
           type - np.array
           desc - score vector
        '''

        idx = np.flatnonzero(sv == np.max(sv))
        if len(idx) > 1:
            print(idx)
            return None, None
        idx = idx[0]
        return idx, sv[idx]

    def predict(self, x):
        '''
        Predict the label of sequence x.

        :: x
           type - string
           desc - sequence to predict
        '''

        # Did we train this instance?
        if not self.primed:
            raise RuntimeError(
                'Unprimed Scry instance')

        # Query length?
        n = len(x)

        # Localize instance attributes
        # for speed (Python is weird!)
        sn = self.n
        st = self.pt

        # Query unresolvable?
        if (n < (sn - st)) or (n > (sn + st)):
            return None, None # Too Many Indels!

        # Case 1: Near Exact Matches

        # Localize Corpus Map
        sc = self.c

        # Quick Check!
        if x in sc:
            y = sc[x]
            if y is None:
                # Uniquely Resolved
                return None, None
            # Pre-resolved Duplicate
            return y, 1.0

        # Analyze Substrings in Decreasing
        # Length Order, Left to Right
        sl = sn if sn < n else n-1
        ys = set()
        pl = sl
        for sx,l in Scry.stream_substrings(
            x=x, p=sl, q=sn-st):
            # Same Length Substrings Analyzed
            if l < pl:
                if len(ys) > 0:
                    if len(ys) == 1:
                        # Uniquely Resolved
                        return ys.pop(), 1.0
                    # Resolution Ambiguous
                    return None, None
                pl = l # Update pl-value
            # Regular Query Analysis
            if sx in sc:
                y = sc[sx]
                if y is None:
                    # Unresolvable Duplicate
                    return None, None
                else:
                    # Potentially Resolved
                    ys.add(y)

        # Last Round Assessment
        if len(ys) > 0:
            if len(ys) == 1:
                # Uniquely Resolved
                return ys.pop(), 1.0
            # Resolution Ambiguous
            return None, None

        # Are we done?
        if self.pm == 0:
            # Continue No Further
            return None, None

        # Case 2: Critical Matches

        # Localize Spectrum Dictionary
        sk = self.k
        sd = self.d
        z  = self.z
        sv = np.zeros(z)
        xs = 0.
        u  = sn-sk

        # Resolve k-mers
        for kmer,idx in Scry.stream_kmers(
            x=x, n=n, k=sk, locations=True):

            # Have we seen this k-mer?
            if kmer in sd:

                # Process Index Range
                iv,ig = sd[kmer]
                mi = idx-st
                si = mi if mi > 0 else 0
                xi = idx+st
                ei = xi if xi < u else u

                # Finalize Index Range
                kk = si
                sr = None
                en = None
                while kk <= ei:
                    if ig[kk] is None:
                        kk += 1
                        continue
                    s,e = ig[kk]
                    if sr is None:
                        sr = s
                    en  = e
                    kk += 1

                # Range Valid?
                if not st is None:
                    # Score Matching Indices
                    xs += sk
                    Scry.update_scorevector(
                        sv=sv, iv=iv, sr=sr, en=en, sk=sk)

        # Did we get any matches?
        if xs == 0.:
            # Nothing Matched ..
            return None, None
        else:
            # Get Best Matches
            idx,idxs = Scry.select_maxindex(
                sv=sv)

            # Ambiguous Matching
            if idx is None or idxs == 0.:
                return None, None

            # True Alignment Scoring
            alns = ed.align(
                query=self.X[idx],
                target=x,
                mode='NW',
                task='distance',
                k=sk)['editDistance']

            # Final Verdict
            if alns == -1:
                return None, None
            return idx, idxs / xs