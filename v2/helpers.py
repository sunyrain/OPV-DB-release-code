from pathlib import Path
import re
import numpy as np
import pandas as pd
EMPTY={'','none','null','nan','n/a','na','unknown','not reported','[]','{}'}
def read(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False, encoding='utf-8-sig')

def present(s):
    return ~s.astype(str).str.strip().str.lower().isin(EMPTY)

def number(s):
    return pd.to_numeric(s, errors='coerce').replace([np.inf,-np.inf],np.nan)

def doi(value):
    return re.sub(r'^https?://(?:dx\.)?doi\.org/','',str(value).strip().lower()).rstrip('.,;')

def label(frame,role):
    return frame[role+'_canonical'].where(present(frame[role+'_canonical']),frame[role]).str.strip()

def family(value):
    if str(value).strip().lower() in EMPTY:
        return 'Unclassified'
    # A reproducible label match, not a structural assignment of every acceptor.
    compact=re.sub(r'[\s\[\](){}_-]','',str(value).upper())
    return 'Fullerene label match' if re.search(r'PCBM|PC6[01]|PC7[01]|ICBA|C60|C70|FULLERENE',compact) else 'Other named acceptors'
