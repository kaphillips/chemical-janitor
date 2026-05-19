import re
import pandas as pd
from .utils import append_col


def re_casrn():
    """
    CASRN Compiled Regular Expression
    """
    return re.compile("[1-9][0-9]{1,6}\-[0-9]{2}\-[0-9]")

def re_formula():
    """
    Molecular Formula Compiled Regular Expression
    """
    s = re.compile('([A-Z][a-z]?)(\d*(?:(?:[\.|\,])\d+(?:\%)?)?)|(?:[\(|\[])([^()]*'
                   '(?:(?:[\(|\[]).*(?:[\)|\]]))?[^()]*)(?:[\)|\]])(\d*(?:(?:[\.|\,]?)'
                   '\d+(?:\%)?))')
    return re.compile(s)


def casrn_split(x):
    """
    Take a text string and extract all occurrences of a CASRN number. This 
    should remove leading zeros (0) too.
    
    Parameters
    ----------
    x: string, text that will be searched for CASRNs
    
    Returns
    -------
    s: integer, a count of how many CASRNs were found in the text
    """
    if isinstance(x,str):
        s = re.findall(re_casrn(),x)
        if len(s) < 1:
            s = pd.NA
    else:
        s = pd.NA
    return s

def check_casrn(x):
    """
    Check that the last digit in the CAS-RN is a valid digit. One way to ensure
    if a CAS-RN is invalid

    iN_i + ... + 4N_4 + 3N_3 + 2N_2 + 1N_1         R
    -------------------------------------- = Q + ----
                       10                         10
                       
    Q is the last digit of the CAS-RN
    
    Parameters
    ----------
    x: string, string of a CAS-RN
    
    Return
    ------
    boolean: True if checksum if valid, False otherwise
    """

    valid = False

    if not isinstance(x,str):
        return valid

    if not re.match(re_casrn(),x):
        return valid

    cas = x[-3::-1].replace('-', '')
    q = 0
    for i,d in enumerate(cas):
        q += (i+1)*int(d)

    if q%10 == int(x[-1]):
        valid = True
    else:
        valid = False
        
    return valid



def find_formula(x):
    """
    Finds a Molecular Formula-like string; this is a stop-gap. It really should
    look for combinations of legitimate elements for a formula, but I haven't
    the time. That's why you see the weird "if" statement to re-include NaCl.
    
    Parameters
    ----------
    x: string in which to search for molecular formulas
    
    Returns
    -------
    boolean, True if there is a formula in the string, False if not, NA if no
    string is passed.
    """
    if isinstance(x,str):
        s = re.findall(re_formula(), x)
        if len(s) < 1:
            s = ''
        else:
            s = "".join(["".join(elem) for elem in s])
            c = re.findall('\d',s)
            if len(c) < 1:
                ## Probably will require manually encoding all of the periodic
                ## table to find all ionic compounds with -1 and +1 charges
                if s != "NaCl":
                    s = ''
    else:
        s = ''
        
    if pd.isnull(x):
        x = 'empty'
    return s == x



def correct_formula(df,col='chemical_name',comment='name_comment'):
    """
    Removes "chemical names" that are only chemical formulas in hiding
    """
    df = df.copy()
    df['name_is_formula'] = df[col].apply(find_formula)
    idx = df.name_is_formula
    df.loc[idx,comment] = (df[idx]
                           .apply(lambda x: append_col(x=x[comment],
                                                       s=x[col],
                                                       comment="Name only formula"),
                               axis=1))
    df.loc[idx,col] = pd.NA
    return df.drop(['name_is_formula'],axis=1).copy()

def string_not_casrn(df,col='casrn',comment='casrn_comment'):
    """
    Remove CAS-RN records that really just contain regular text
    """
    df = df.copy()
    idx = ~df[col].str.contains("[1-9][0-9]{1,6}\-[0-9]{2}\-[0-9]",regex=True,na=True)
    df.loc[idx,comment] = (df[idx]
                                .apply(lambda x: append_col(x[comment],
                                                            s=x[col],
                                                            comment="String is not CAS-RN"),
                                        axis=1))
    df.loc[idx,col] = pd.NA
    return df.copy()



def split_casrns(df,col='casrn',comment='casrn_comment'):
    """
    If multiple CAS-RNs are on one line, then split them
    """
    df = df.copy()
    df[col] = df[col].apply(casrn_split)
    df[f'{col}_len'] = df[col].str.len()
    idx = (df[f'{col}_len'] != 1) & (df[f'{col}_len'].notnull())

    df.loc[idx,comment] = (df[idx]
                        .apply(lambda x: append_col(x[comment],
                                                    s=",".join(x[col]),
                                                    comment="Multiple CAS-RN on line"),
                                axis=1))
    df.drop([f'{col}_len'],axis=1,inplace=True)
    df = df.explode(column=col).copy()
    return df.copy()



def casrn_checksum(df,col='casrn',comment='casrn_comment'):
    """
    Perform CAS-RN checksum
    """
    df = df.copy()
    idx = (~df[col].apply(check_casrn)) & (df[col].notnull())
    df.loc[idx,comment] = (df[idx]
                        .apply(lambda x: append_col(x[comment],
                                                    s=x[col],
                                                    comment="CAS-RN failed checksum"),
                                axis=1))
    df.loc[idx,col] = pd.NA
    return df.copy()


def casrn_finder(df,cas_col='casrn',name_col='chemical_name',comment='name_comment'):
    
    df = df.copy()
    df['cas_in_name'] = df[name_col].apply(casrn_split)
    
    idx = df['cas_in_name'].notnull()
    
    df.loc[idx,'cas_in_name'] = df.loc[idx,'cas_in_name'].apply(lambda x: ", ".join(x))
    
    df[cas_col] = (df[[cas_col,'cas_in_name']]
                   .apply(lambda x: ", ".join(set([i for i in x if pd.notnull(i)])), 
                          axis=1))

    df.loc[idx, comment] = (df[idx]
                            .apply(lambda x: append_col(x=x[comment],
                                                        s=x[cas_col],
                                                        comment="CAS-RN in name; copied to casrn column"),
                                   axis=1))

    regex = re.compile("\(CAS Reg. No. [1-9][0-9]{1,6}\-[0-9]{2}\-[0-9]\)")
    df.loc[idx,name_col] = (df.loc[idx,name_col]
                            .apply(lambda x: re.sub(regex,"",x)))
    
    return df.copy()