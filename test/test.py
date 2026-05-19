import pandas as pd
from chemical_janitor.chemicals import (
    correct_formula,
    string_not_casrn,
    split_casrns,
    casrn_checksum,
    casrn_finder,
)
from chemical_janitor.drops import(
    drop_terminal_phrases,
    drop_fcs,
    drop_stoppers,
    drop_foods,
    drop_text,
    drop_salts,
    drop_blocks
    )

from chemical_janitor.encodings import fix_encodings
from chemical_janitor.termini import terminal_unspecified
from chemical_janitor.utils import string_cleaning, newest_file, date_file

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

