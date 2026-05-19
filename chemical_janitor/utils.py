import string
from pathlib import Path
from datetime import datetime

import pandas as pd



def append_col(x,s,comment,sep="|"):
    """
    Used for appending reasons why chemical name or casrn was altered. 
    
    Parameter
    ---------
    x: None or string; original comment
    s: None or string; information to append to original comment
    comment: string; context to add to the new information
    sep: character to separate the various comments
    
    Returns
    -------
    None if there is no informaiton to add; String of new information +
    old information, otherwise
    """
    if isinstance(x,str):
        if isinstance(s,str):
            s = f"{comment}: {s}"
            y = sep.join([x.strip(),s.strip()])
        else:
            y = x
    elif pd.isnull(x):
        if pd.notnull(s):
            s = f"{comment}: {s.strip()}"
        y = s
    else:
        y = pd.NA
        
    return y



def newest_file(globber,path=Path()):
    """
    Returns the newest file that matches a glob statement based on the "creation
    time"
    path: string, Path; defaults to .: directory in which to look for file
    globber: string: wild-card search string to use for finding matching files
    """
    if isinstance(path,str):
        path = Path(path)
    path = path.glob(globber)
    return max(path,key=lambda x: x.stat().st_ctime)



def date_file(stem, suffix, sep="-", format='%b-%d-%Y'):
    """
    Introduce a date stamp into a string.

    Parameters
    ----------
    stem : the stem of file
    suffix : the extension of the file
    out : a string that is of the form prefix_MMDDYYYY.suffix
    """
    hoy = datetime.date.today().strftime(format)
    
    return f"{stem}{sep}{hoy}.{suffix.strip('.')}"


def string_cleaning(df,col):
    """
    Remove leading or trailing white spaces or punctuation
    """
    df = df.copy()
    omits = (string.whitespace+
             (string.punctuation
              .replace("()","")
              .replace("[","")
              .replace("]","")
              .replace("{","")
              .replace("}",""))+
             string.whitespace)

    for p in omits:
        df[col] = df[col].str.strip(p)
    df[col] = df[col].str.strip()
    
    return df.copy()



