import numpy as np
import pandas as pd

from .chemicals import casrn_split


from pandas.api.extensions import register_series_accessor

def _changed(before: pd.Series, after: pd.Series) -> pd.Series:
    return ~before.fillna("").eq(after.fillna(""))


def _append_note(notes: pd.Series, mask: pd.Series, msg: str) -> pd.Series:
    notes = notes.astype("object").fillna("")
    addition = pd.Series(np.where(mask, msg, ""), index=notes.index, dtype="object")
    out = notes.where(
        addition.eq(""),
        notes.str.cat(addition, sep="; ").str.strip("; ")
    )
    return out.replace("", pd.NA)

@register_series_accessor("casrn")
class CASRNAccessor:
    """
    CAS-RN Validation and Cleaning
    """
    
    def __init__(self, obj: pd.Series):
        self._obj = obj.astype(pd.StringDtype())
        self._notes = pd.Series(index=obj.index,dtype=pd.StringDtype())

    def _note(self, mask: pd.Series,  msg: str):
        self._notes = _append_note(self._notes,mask,msg)

    def not_casrn(self):
        before = self._obj
        idx = ~before.str.contains("[1-9][0-9]{1,6}\-[0-9]{2}\-[0-9]",regex=True,na=True)
        after = before.copy()
        after[idx] = pd.NA
        self._note(_changed(before,after),"String is not CAS-RN")
        self._obj = after
        return self

    def split(self):
        before = self._obj
        after = before.apply(casrn_split)
        casrn_len = after.str.len()
        idx = (casrn_len != 1) & (casrn_len.notnull())
        





