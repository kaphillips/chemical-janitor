import re

import pandas as pd

from .stoppers import function_categories, fc_string, foods, standard_stops, block_list
from .termini import terminal_bracket, terminal_parenthesis
from .utils import append_col

## FIXME: how does block_list differ from standard_stops? I should make a class that
## has different types of stop words (functional use, foods, standard stops, and block
## list, custom list) and allow the user the to select the type of stop words they'd
## like implemented in their workflow


def drop_fcs(df,col='chemical_name',comment='name_comment'):
    """
    Drop chemical names that are just FCs
    """
    df = df.copy()
    
    ## All the cleaned function categories
    fcs = function_categories().fc_clean
    idx = df[col].str.lower().isin(fcs)
    df.loc[idx,comment] = (df[idx]
                               .apply(lambda x: append_col(x=x[comment],
                                                           s=x[col],
                                               comment="Name is functional use"),
                           axis=1))
    df.loc[idx,col] = pd.NA

    idx = ((df[col].str.lower().str.contains(fc_string(),regex=True,na=False)) &
        (~df[col].str.contains("\d",regex=True,na=False)))
    
    df.loc[idx,comment] = (df[idx]
                           .apply(lambda x: append_col(x=x[comment],
                                                       s=x[col],
                                           comment="Name is functional use"),
                           axis=1))
    df.loc[idx,col] = pd.NA
    return df.copy()



def drop_foods(df,col='chemical_name',comment='name_comment'):
    """
    Drop chemical names that are just foods
    """
    df = df.copy()
    idx = (df[col].str.lower().str.contains(foods(),regex=True,na=False))
    df.loc[idx,comment] = (df[idx]
                           .apply(lambda x: append_col(x=x[comment],
                                                       s=x[col],
                                              comment="Name is food"),
                           axis=1))
    df.loc[idx,col] = pd.NA
    return df.copy()



def drop_stoppers(df,col='chemical_name',comment='name_comment'):
    """
    Drop stop words, 
    """
    df = df.copy()

    idx = ((df[col]
            .str.lower()
            .str.contains("|".join(standard_stops()),
                          regex=True,na=False)) & 
           (~df[col]
            .str.lower()
            .str.contains("yl",
                          regex=False,na=False)))
    df.loc[idx,comment] = (df[idx]
                                    .apply(lambda x: append_col(x=x[comment],
                                                                s=x[col],
                                                       comment="Ambiguous name"),
                                    axis=1))
    df.loc[idx,col] = pd.NA

    idx = df[col].str.lower().isin(["polymer",'polymers','wax',"mixture"])
    df.loc[idx,comment] = (df[idx]
                           .apply(lambda x: append_col(x=x[comment],
                                                       s=x[col],
                                              comment="Ambiguous name"),
                           axis=1))
    df.loc[idx,col] = pd.NA

    ## The stops() function gets rid of some, but the "yl" keeps a few, remove them
    ## here
    idx = df[col].str.lower().str.contains("citron",na=False)
    df.loc[idx,comment] = (df[idx]
                           .apply(lambda x: append_col(x=x[comment],
                                                       s=x[col],
                                              comment="Ambiguous name"),
                           axis=1))
    df.loc[idx,col] = pd.NA
    
    idx = df[col].str.lower().str.contains("compound",na=False)
    df.loc[idx,comment] = (df[idx]
                           .apply(lambda x: append_col(x=x[comment],
                                                       s=x[col],
                                              comment="Ambiguous name"),
                           axis=1))
    df.loc[idx,col] = pd.NA
    return df.copy()



def drop_text(df,col='chemical_name',comment='name_comment'):
    """
    Drop extraneous text (various phrases that occurred frequently)
    """
    df = df.copy()
    
    ## If a string starts with "Part *:"
    idx = (df[col].str.lower().str.contains("part [a-z]:",regex=True,na=False))
    df.loc[idx,col] = df.loc[idx,col].str.split(":")
    df.loc[idx,comment] = (df[idx]
                           .apply(lambda x: append_col(x=x[comment],
                                                       s=x[col][0],
                                           comment="Removed text"),
                           axis=1))
    df.loc[idx,col] = df.loc[idx,col].str[-1]
    
    ## If a chemical name has "modified" in some form in the name, remove it
    idx = df[col].str.lower().str.contains('modif',na=False,regex=True)
    df.loc[idx,comment] = (df[idx]
                           .apply(lambda x: append_col(x=x[comment],
                                                       s=x[col],
                                           comment="Unknown modification"),
                           axis=1))
    df.loc[idx,col] = pd.NA

    ## If a chemical name has "pure", 
    quality = ['pure','purif','tech','grade','chemical']
    pat = [fr'(\w*{word}\w*)' for word in quality]
    pat = fr'{"|".join(pat)}'
    pat = re.compile(pattern=pat,flags=re.IGNORECASE)
    idx = df[col].str.lower().str.contains('|'.join(quality),na=False,regex=True)
    df.loc[idx,comment] = (df[idx]
                           .apply(lambda x: append_col(x=x[comment],
                                                       s=" ".join([j for tup in re.findall(pattern=pat,
                                                                              string=x[col]) for j in tup if j != ""]),
                                           comment="Unneeded adjective"),
                           axis=1))

    df.loc[idx,col] = (df.loc[idx,col].apply(lambda x: re.sub(pattern=pat,
                                                                repl="",
                                                                string=x)))

    ## Terminal percentage
    pat = re.compile('\d+\%$',flags=re.IGNORECASE)
    idx = df[col].str.lower().str.contains('.*\d\%$',na=False,regex=True)
    df.loc[idx,comment] = (df[idx]
                           .apply(lambda x: append_col(x=x[comment],
                                                       s=(re.search(pattern=pat,
                                                                   string=x[col])
                                                          .group()
                                                          .strip()),
                                           comment="Removed text"),
                           axis=1))
    df.loc[idx,col] = (df.loc[idx,col].apply(lambda x: re.sub(pattern=pat,
                                                                repl="",
                                                                string=x)))
    df[col] = df[col].str.strip().str.strip(',').str.strip('-').str.strip()
    df[comment] = df[comment].str.strip()
    return df.copy()



def drop_salts(df,col='chemical_name',comment='name_comment'):
    """
    Drop reference to ambiguous salts
    """
    df = df.copy()
    pat = re.compile('and its .* salts|and its salts',flags=re.IGNORECASE)
    idx = df[col].str.lower().str.contains('and its .* salts|and its salts',regex=True,na=False)
    df.loc[idx,comment] = (df[idx]
                            .apply(lambda x: append_col(x=x[comment],
                                                        s=(re.search(pattern=pat,
                                                                     string=x[col])
                                                            .group()
                                                            .strip()),
                                            comment="Ambiguous salt reference"),
                            axis=1))
    df.loc[idx,col] = (df.loc[idx,col]
                       .apply(lambda x: re.split(pattern=pat,
                                                string=x,)[0]))
    return df.copy()




def drop_blocks(df,col='chemical_name',comment="name_comment"):
    """
    Drop block words
    """
    df = df.copy()
    blocks = block_list()
    
    idx = (df[col].isin(blocks))
    
    df.loc[idx,comment] = (df[idx]
                           .apply(lambda x: append_col(x=x[comment],
                                                       s=x[col],
                                                       comment="Name is on block list"),
                                  axis=1))
    
    df.loc[idx,col] = pd.NA

    return df.copy()



def drop_terminal_phrases(df, col='chemical_name', comment='name_comment'):
    """
    Drop terminal parenthesis or brackets
    """
    
    df = df.copy()
    df['parenth'] = df[col].apply(terminal_parenthesis)
    df[comment] = (df
                   .apply(lambda x: append_col(x=x[comment],
                                               s=x['parenth'][-1],
                                               comment="Extraneous parenthesis"),
                          axis=1))
    df[col] = df['parenth'].str[0]

    df['brackets'] = df.chemical_name.apply(terminal_bracket)
    df[comment] = (df
                   .apply(lambda x: append_col(x=x[comment],
                                               s=x['brackets'][-1],
                                               comment="Extraneous brackets"),
                          axis=1))
    df[col] = df['brackets'].str[0]
    df[col] = df[col].str.strip()

    return df.drop(['parenth', 'brackets'], axis=1).copy()