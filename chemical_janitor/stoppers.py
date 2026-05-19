import pandas as pd
from ctxpy import Exposure


def standard_stops():
    """
    Stop-words that have been compiled based on known reported "chemical names" that
    are not really chemical names

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


def custom_stops(path):
    """
    Read in a custom list of stop-words from a text file: each line in a text file is 
    a single stop word. Regular Expression characters are accepted.
    """
    with open(path,'r') as f:
        stop_words = f.readlines()
    return stop_words



def function_categories():
    """
    Don't keep any chemical names whose "name" matches a functional use category
    from the OECD/EPA
    """

    # Get current list of CPDat Function Categories
    df = (Exposure()
          .vocabulary(by='fc').to_df()
          .drop(['id','description'],axis=1)
          .rename(columns={"title":"function_category"}))

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