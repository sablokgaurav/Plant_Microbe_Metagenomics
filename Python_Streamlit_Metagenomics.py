# a streamlit for searching the metagenomics abundance species
# for training a language model for the abundant species 
# final check  
import streamlit as st
import pandas as pd
import wikipedia as wikipedia
import csv 
from googlesearch.googlesearch import GoogleSearch
st.set_page_config(

    page_title="Welcome to the metagenomics miner",
    page_icon="I am going to search for your species",
    layout="wide",
    initial_sidebar_state="expanded",

)
st.title("Metagenomics Dataframe Processing")
st.header('''You can search the species from your abundance 
                classification or you can add your 
                               own species to the classified datasets''')
datframe = st.file_uploader("Please upload your metagenomics dataframe:")
dataframeU = st.dataframe(dataframe)
speciesCol = st.text("Enter the species column name:")
st.button ("Click for the processing of the dataframe")
if st.button:
    selected_species = list(set(dataframeU["speciesCol"].dropna().tolist()))
st.text("Metagenomics search engine")
search_engine = st.multiple("Select your search engine", ["Google", "Wikipedia"])
st.button("Choose the search engine")
if st.button and search_engine == "Google":
    selected_species = list(set(dataframeU["speciesCol"].dropna().tolist()))
    for i in range(len(selected_species)):
        species_text = GoogleSearch().search(selected_species[i])
        st.write(f'''Your bacterial species title is:
                      {selected_species[i].title}'''.format)
        st.write(f'''Your selected bacterial species text is: 
                             {selected_species[i].gettext()}''')
if st.button and search_engine == "Wikipedia":
    for i in range(len(selected_species)):
        st.write(f'''the wikipedia species summary is: 
                      {wikipedia.summary(selected_species[i]), sentences == int(101)}''')
        st.write(f'''the wikipedia species search is
                                        :{wikipedia.search(selected_species[i])}''')
        st.write(f'''the_wikipedia page for the species is 
                                         :{wikipedia.page(selected_species[i])}''')
else:
    print("No species selected")
    
additional_species = st.text_area("Add additional species to the existing datasets")
st.button("Adding additional species to the existing search")
if st.button:
    addition = list(set(additional_species.strip().split()))
    total = selected_species+addition
    st.write(f'the total number of the species for the search are:{total}')
st.button("Search entire species")
if st.button and search_engine == "Google":
     total_species = selected_species+addition
for i in range(len(total_species)):
        species_text = GoogleSearch().search(total_species[i])
        st.write(f'Your bacterial species title is:{total_species[i].title}')
        st.write(f'''Your selected bacterial species text is: 
                                             {total_species[i].gettext()}''')
if st.button and search_engine == "Wikipedia":
        for i in range(len(total_species)):
            st.write(f'''the wikipedia species summary is: 
                        {wikipedia.summary(total_species[i]), sentences = 101}''')
            st.write(f'''the wikipedia species search is: 
                                          {wikipedia.search(total_species[i])}''')
            t.write(f'''the_wikipedia page for the species is 
                                          :{wikipedia.page(total_species[i])}''')
else:
    print("No search selected")
   
st.text("I can help you also perform a selected search from a dataframe")
selected_species = st.select_slider("select a species from a given datafram",
                    list(set(dataframeU["speciesCol"].dropna().tolist())))
engine = st.multiple("Select your search engine": ["Google", "Wikipedia"])
st.button("Search selected species")
if st.button and engine=="Google":
    species_text = GoogleSearch().search(selected_species)
    st.write(f'the title of the species is: {species_text.title}')
    st.write(f'the text of the species is: {species_text.gettext()} ')
elif st.button and engine == "Wikipedia":
    st.write(f'''the wikipedia species summary is: 
                   {wikipedia.summary(total_species[i]), sentences == 101}''')
    st.write(f'''the wikipedia species search is: 
                      {wikipedia.search(total_species[i])}''')
    st.write(f'''the_wikipedia page for the species is 
                                      :{wikipedia.page(total_species[i])}''')
else:
    print("No search selected")
