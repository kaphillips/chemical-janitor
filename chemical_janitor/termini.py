import re
import pandas as pd
from .utils import append_col

def terminal_parenthesis(x):
    """
    Remove terminal parenthesis if it doesn't contain part of a chemical name
    """
    if isinstance(x,str):
        ## Returns only last match...
        s = re.findall('.*\((.*?)\)',x)
        if len(s) < 1:
            s = pd.NA
        else:
            ## ...but it's still a list so, pull it out
            s = s[-1]
            
            ## If the text in the parenthesis isn't at the end of the string,
            ## don't remove it, exit search
            if x[-(len(s)+1):-1] != s:
                s = pd.NA
            
            ## A lot of chemicals have "yl" in the string, yet it is not a
            ## common letter combination seen in the rest of the English 
            ## language use this to find as many last parenthetical phrases 
            ## that contain a chemical name (and therefore shouldn't be removed)
            ## as many as possible
            if pd.notnull(s):
                if "yl" in s:
                    keepers = ['density','probably','average','combination']
                    if not any(i in keepers for i in s.split()):
                        s = pd.NA
                
    else:
        s = pd.NA
        

    if pd.isnull(s):
        phrase = x
    else:
        phrase = x[:-(len(s)+2)].strip()
    return (phrase,s)



def terminal_bracket(x):
    """
    Remove terminal bracket if it doesn't contain part of a chemical name
    """
    if isinstance(x,str):
        ## Returns only last match...
        s = re.findall('.*\[(.*?)\]',x)
        if len(s) < 1:
            s = pd.NA
        else:
            ## ...but it's still a list so, pull it out
            s = s[-1]
            
            ## If the text in the parenthesis isn't at the end of the string,
            ## don't remove it, exit search
            if x[-(len(s)+1):-1] != s:
                s = pd.NA
            
            ## A lot of chemicals have "yl" in the string, yet it is not a
            ## common letter combination seen in the rest of the English 
            ## language use this to find as many last parenthetical phrases 
            ## that contain a chemical name (and therefore shouldn't be removed)
            ## as many as possible
            if pd.notnull(s):
                if "yl" in s:
                    keepers = ['density','probably','average','combination']
                    if not any(i in keepers for i in s.split()):
                        s = pd.NA
                
    else:
        s = pd.NA
        

    if pd.isnull(s):
        phrase = x
    else:
        phrase = x[:-(len(s)+2)].strip()
    return (phrase,s)



def terminal_unspecified(df,col='chemical_name',comment='name_comment'):
    """
    Removes terminal "{PUNCT} unspecified" phrase from the end of strings
    """
    df = df.copy()
    idx = (df[col]
           .str.lower()
           .str.contains('[.?\-",]+ unspecified',
                         na=False,regex=True))
    
    df.loc[idx,comment] = (df[idx]
                           .apply(lambda x: append_col(x=x[comment],
                                                       s=re.search('[.?\-",]+ unspecified',
                                                                   x[col],
                                                                   flags=re.IGNORECASE).group().strip(),
                                           comment="Unspecified warning"),
                           axis=1))
    df.loc[idx,col] = (df
                       .loc[idx,col]
                       .apply(lambda x: re.sub('[.?\-",]+ unspecified',"",
                                               x.strip(",").strip(),
                                               flags=re.IGNORECASE)))
    return df.copy()
