from .utils import append_col

def has_unicode(x):
    """
    Checks if there is a unicode character (or character) in the passed string
    
    Parameters
    ----------
    x: string to check for unicode
    
    Returns
    -------
    boolean, True if there is at least 1 unicode character; False if not
    """
    if isinstance(x,str):
        enc = x.encode('utf-8',errors='replace')
        dec = enc.decode('utf-8')
        return len(enc)!=len(dec)
    else:
        return False



def known_encodings():
    """
    UTF-8 code points for which I can sub in ASCII text
    """
    ## Can be added to, if needed
    unis = {"\u2032": "'",
            "\u03c9": ".omega.",
            "\xae": " (registered trademark)",
            "\u2013": "--",
            "\xb0": " degrees ",
            "\u2019": "'",
            "\u2026": "...",
            "\u03b1": ".alpha.", }
    return unis

def fix_encodings(df,col='chemical_name',comment='name_comment'):
    """
    Swap out known UTF-8 encodings with ASCII
    """
    df = df.copy()

    # Dictionary of unicode codepoints and ascii replacements
    unis = known_encodings()
    
    # Find which rows in name and cas have unicode characters
    df['has_unicode'] = df[col].apply(has_unicode)

    ## Loop over each known unicode codepoint
    for k,v in unis.items():
        
        idx = df[col].str.contains(k,na=False)
        df[col] = df[col].str.replace(k,v)
        df.loc[idx,comment] = (df[idx]
                               .apply(lambda x: append_col(x[comment],
                                                           s=f"swapped {k} with {v}",
                                                           comment="Unicode detected"),
                                      axis=1))

    df.drop(['has_unicode'],axis=1,inplace=True)
    return df.copy()