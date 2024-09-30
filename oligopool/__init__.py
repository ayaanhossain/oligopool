# Core Design Functions
from .base import background
from .base import barcode
from .base import primer
from .base import motif
from .base import spacer

# Assembly Design Functions
from .base import split
from .base import pad

# Auxiliary Design Functions
from .base import lenstat
from .base import final

# Analysis Functions
from .base import index
from .base import pack
from .base import acount
from .base import xcount

# Utilities
from .base import utils

# Aliasing
background = background.background
barcode = barcode.barcode
primer = primer.primer
motif = motif.motif
spacer = spacer.spacer
split = split.split
pad = pad.pad
lenstat = lenstat.lenstat
final = final.final
index = index.index
pack = pack.pack
acount = acount.acount
xcount = xcount.xcount

# Setup
__author__ = 'Ayaan Hossain'

__version__ = '0.0.0'

__all__ = [
    background, barcode, primer, motif, spacer,
    split, pad, lenstat, final,
    index, pack, acount, xcount,
    utils
]

__doc__ = 'TBW'