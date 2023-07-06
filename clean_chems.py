import string
import pandas as pd
import re
import datetime
from pathlib import Path



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
        s = re.findall("[1-9][0-9]{1,6}\-[0-9]{2}\-[0-9]",x)
        if len(s) < 1:
            s = pd.NA
    else:
        s = pd.NA
    return s



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
        regex = re.compile('([A-Z][a-z]?)(\d*(?:(?:[\.|\,])\d+(?:\%)?)?)|(?:[\(|\[])([^()]*(?:(?:[\(|\[]).*(?:[\)|\]]))?[^()]*)(?:[\)|\]])(\d*(?:(?:[\.|\,]?)\d+(?:\%)?))')
        s = re.findall(regex, x)
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



def checksum(x):
    """
    Check that the last digit in the CAS-RN is a valid digit. One way to ensure
    if a CAS-RN is invalid
    
    Parameters
    ----------
    x: string, string of a CAS-RN
    
    Return
    ------
    boolean: True if checksum if valid, False otherwise
    """
    if not isinstance(x,str):
        return False
    cas = x[-3::-1].replace('-', '')
    q = 0
    for i,d in enumerate(cas):
        q += (i+1)*int(d)
    if q%10 == int(x[-1]):
        return True
    else:
        return False



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



def stops():
    """
    Custom stop-words for automatically removing a "chemical name"

    Returns
    -------
    List of stop-word strings
    """
    stop_words = ['proprietary','ingredient','hazard','blend','inert','stain'
                  'other', 'withheld', 'cas |cas-|casrn', 'secret', "herbal",
                  'confidential','bacteri','treatment','contracept','emission',
                  "agent","eye","resin","citron",'bio','smoke','fiber','adult',
                  'boy','girl','infant','child','other organosilane','material']
    return stop_words



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



def term_parenth(x):
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



def term_bracket(x):
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



def function_categories():
    """
    Don't keep any chemical names whose "name" matches a functional use category
    from the OECD/EPA
    """
    # If chemical name is a functional use category, it's not a real chemical name
    df = (pd.read_csv("ChemExpo_FC_2023-05-24.csv", usecols=[0])
          .rename(lambda x: x.lower().replace(" ", "_"), axis=1))

    # Some FC have terminal parentheticals. Some of these are FC-synonyms, some
    # just denote the FC as an EPA-added FC. Remove the EPA-added parenthetical
    # and join the synonym in a regex-friendly way to facilitate a str.contains
    # statement that looks of any FC (or synonym) in a string
    df['fc_clean'] = (df.function_category
                      .str.lower()
                      .str.split("("))
    df['fc_extra'] = df.fc_clean.str[-1].str.replace(")", "", regex=False)
    df['fc_clean'] = df.fc_clean.str[0]
    df.loc[df.fc_clean == df.fc_extra, 'fc_extra'] = pd.NA
    df.loc[df.fc_extra == "epa", 'fc_extra'] = pd.NA
    df['fc_clean'] = df[['fc_clean', 'fc_extra']].apply(
        lambda x: "|".join([i.strip() for i in x if pd.notnull(i)]), axis=1)
    df.drop(['fc_extra'], axis=1, inplace=True)
    
    return df.copy()



def fc_string():
    """
    This helps account for those FC synoyms
    """
    ## Make a regex string to search for a string containing any of the FCs
    extra_fcs = ['colorant','detergent','additive','flavor',"anti-",
                 "protectant",'thermoplastic',"dispersion","plast",
                 'enzyme','thick','inhib']
    oecd_fcs = function_categories()
    
    extra_fcs = "|".join(extra_fcs)
    oecd_fcs = "|".join(oecd_fcs.fc_clean.unique())
    
    return "|".join([oecd_fcs,extra_fcs])



def foods():
    """
    Foods aren't chemical names either
    """
    food = ['yeast culture', 'food starch', 'sweet whey',
            'salted fish','beverage']
    return "|".join(food)



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


def drop_terminal_phrases(df, col='chemical_name', comment='name_comment'):
    """
    Drop terminal parenthesis or brackets
    """
    
    df = df.copy()
    df['parenth'] = df[col].apply(term_parenth)
    df[comment] = (df
                   .apply(lambda x: append_col(x=x[comment],
                                               s=x['parenth'][-1],
                                               comment="Extraneous parenthesis"),
                          axis=1))
    df[col] = df['parenth'].str[0]

    df['brackets'] = df.chemical_name.apply(term_bracket)
    df[comment] = (df
                   .apply(lambda x: append_col(x=x[comment],
                                               s=x['brackets'][-1],
                                               comment="Extraneous brackets"),
                          axis=1))
    df[col] = df['brackets'].str[0]
    df[col] = df[col].str.strip()

    return df.drop(['parenth', 'brackets'], axis=1).copy()



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
            .str.contains("|".join(stops()),
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
    idx = data[col].str.lower().str.contains('modif',na=False,regex=True)
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
    idx = data[col].str.lower().str.contains('.*\d\%$',na=False,regex=True)
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
    df = data.copy()

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
    idx = (~df[col].apply(checksum)) & (df[col].notnull())
    df.loc[idx,comment] = (df[idx]
                        .apply(lambda x: append_col(x[comment],
                                                    s=x[col],
                                                    comment="CAS-RN failed checksum"),
                                axis=1))
    df.loc[idx,col] = pd.NA
    return df.copy()


def block_list():
    block = ['alcohol', 'alcohol', 'Bly', 'Bly', 'Polyester', 'Polyester',
             'Alkanes', 'Alkanes', 'alkanes', 'alkanes', 'red 4, 33',
             'red 4, 33', 'rose', 'rose',
             'Organic electrolyte principally involves ester carbonate',
             'Organic electrolyte principally involves ester carbonate',
             'PP', 'PP', 'Amine soap', 'Amine soap', 'Free Amines',
             'Free Amines', 'Acrylic Polymer', 'Acrylic Polymers',
             'Urethane Polymer', 'Acrylic Polymer', 'Acrylic Polymers',
             'Urethane Polymer', 'Caustic Salt', 'Caustic Salt', '','',
             'Aflatoxins', 'Aflatoxins', 'Aminoglycosides', 'Anabolic steroids',
             'Analgesic mixtures containing Phenacetin', 'Aminoglycosides',
             'Anabolic steroids', 'Analgesic mixtures containing Phenacetin',
             'Aristolochic acids', 'Aristolochic acids', 'Barbiturates',
             'Barbiturates', 'Benzodiazepines', 'Benzodiazepines',
             'Conjugated estrogens', 'Conjugated estrogens',
             'Dibenzanthracenes', 'Dibenzanthracenes', 'Estrogens, steroidal',
             'Estrogen-progestogen (combined) used as menopausal therapy',
             'Estrogens, steroidal',
             'Estrogen-progestogen (combined) used as menopausal therapy',
             'Etoposide in combination with cisplatin and bleomycin',
             'Etoposide in combination with cisplatin and bleomycin',
             'Cyanide salts that readily dissociate in solution (expressed as cyanide)f',
             'Cyanide salts that readily dissociate in solution (expressed as cyanide)f',]
    return list(set(block))


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



# Load file
ifile = newest_file(globber="uncurated_chemicals*.csv")
data = pd.read_csv(ifile)
data['chemical_name'] = data['raw_chem_name'].copy().astype(pd.StringDtype())
data['casrn'] = data.raw_cas.copy().astype(pd.StringDtype())
data['casrn_comment'] = pd.NA
data['name_comment'] = pd.NA

## String canonicalization for CASRNs
data['casrn'] = (data.casrn
                       .str.strip()
                       .str.replace("...","",regex=False)
                       .str.replace(" (registered trademark)","",regex=False)
                       .str.replace("#","",regex=False)
                       .str.strip("*")
                       .str.split()
                       .str.join(' ')
                       .str.replace(" ",""))

## String canonicalization for Chemical Names
data['chemical_name'] = (data.chemical_name
                           .str.strip()
                           .str.replace("...","",regex=False)
                           .str.replace(" (registered trademark)","",regex=False)
                           .str.replace("#","",regex=False)
                           .str.strip("*")
                           .str.split()
                           .str.join(' '))


## Swap empty CASRNs for Nulls
data.loc[data.casrn=='-','casrn'] = pd.NA

## Chemical Name cleaning steps
data = fix_encodings(df=data)
data = correct_formula(df=data)
data = drop_terminal_phrases(df=data)
data = drop_fcs(df=data)
data = drop_foods(df=data)
data = drop_stoppers(df=data)
data = drop_text(df=data)
data = terminal_unspecified(df=data)
data = drop_salts(df=data)
data = drop_terminal_phrases(df=data)
data = string_cleaning(df=data,col='chemical_name')

## Find CASRNs in chemical name and move to CASRN column
data = casrn_finder(df=data)


## CASRN Cleaning Steps
data = string_not_casrn(df=data)
data = split_casrns(df=data)
data = casrn_checksum(df=data)
data = string_cleaning(df=data,col='casrn')
data['chemical_name'] = data['chemical_name'].str.split().str.join(" ")
data = drop_blocks(df=data)

## Drop records where both clean name and clean casrn is null
all_null = (data.chemical_name.isnull()) & (data.casrn.isnull())
data = data[~all_null].copy()

## This is a manual flag, on 6/7/23 AJW pointing out that there are records 
## whose chemical names have "cyandidef". Saskhi looked it up on 6/30/2023 and
## confirmed that this was an extraction error and that these records should be
## removed. We should look into this more in the future.
data = data[~data.chemical_name.str.contains("cyanidef",na=False)].copy()


data.to_excel(date_file("cleaned_chemicals_for_curation","xlsx"),index=False)





##TODO
## carbon backbone lengths ex. C8-C10
## Greek letters: are they preferred to be .alpha.
## What to do when CASRNs are in name?
## "with the exception of those specified elsewhere in this annex"
## used as; combination, mix
## metabolize, transform, react, salts, by-product
## colon separators
## components
## aromatics, percentages
## Strip punctuation at the end (and %), plus white space, plus reduce double white space
