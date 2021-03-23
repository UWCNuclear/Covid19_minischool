# prerequisite packages to install
# python3
# pip3 install requests
# pip3 install numpy
# pip3 install pandas
# pip3 install setuptools
# pip3 install bs4
# pip3 install tqdm

#####################################
## include needed python utilities ##
#####################################

import requests
import re
import os
import numpy as np
import glob
import pandas as pd
import time

from bs4 import BeautifulSoup
from tqdm import tqdm
from time import gmtime, strftime

###########
## Paths ##
###########

# where the data will be stored
datapath = "./worldometers"

# The web site of worldometers
baseurl = "https://www.worldometers.info/coronavirus/"

#########################################################################
## Get the total list of countries available that we want to update    ##
## if selected_countries is set to "", all countries will be downloaed ##
## => All means 218 countries, can be long                             ##
#########################################################################

#selected_countries = "USA ; Brazil ; UK ; Mexico ; Italy ; France ; Spain ; India ; Iran ; Peru ; Russia ; Belgium ; Germany ; Canada ; Chile ; Netherlands ; Colombia ; Sweden ; Turkey ; Pakistan ; Ecuador ; South Africa ; China"
selected_countries = "USA ; Brazil ; South Africa"
#selected_countries=""

#######################################
## main, where the code is executed  ##
#######################################

def main():

    # we start the time counter
    start = time.time()

    # create, if don't exists, the output directory
    global datapath
    if not datapath.endswith("/"):
        datapath += "/"
    if not os.path.exists(datapath):
        os.makedirs(datapath)
    print("The files are downloaded in: " + datapath)


    # if not existing, we build the file containing all the available countries
    if not os.path.isfile(os.path.join(datapath, 'country_list.csv')):
        df_country_list = get_country_list()
        df_country_list.to_csv(os.path.join(datapath, 'country_list.csv'))

    # we read the file containing all the countries, and each associated link
    df_country_list = pd.read_csv(os.path.join(datapath, 'country_list.csv'))[['country','link']]

    # now, we loop on all the requested countries
    for i, row in df_country_list.iterrows():

        c_name = row['country']
        link = row['link']

        # if all the countries are requested, or if the current country is in the requested list, we process the download, else, we go the next country
        if(selected_countries == "" or selected_countries.find(c_name) >= 0):
            print(i, c_name, end=" ")
        else:
            continue

        # we create a DataFrame: Two-dimensional, size-mutable, potentially heterogeneous tabular data, to store the data
        df = pd.DataFrame()

        # we parse the content of the htmp web page of the specific country
        link = baseurl+link
        try:
            html_page = load_data(link)
            res = html_parser(html_page.content, 'script')
        except Exception as e:
            print(" -> oups, it seems to have an issue with that country :( :( :( ")
            continue
        print(".",end="")

        # then, we extract the dates and data contained in the graph of Total Cases
        date_range, data_total_cases = clean_data(res, 'Total Cases')
        if date_range != False:
            df['Date'] = date_range
            df['Total Cases'] = data_total_cases

        # same for Total Deaths
        date_range, data_total_deaths = clean_data(res, 'Total Deaths')
        if date_range != False:
            df['Total Deaths'] = data_total_deaths

        # we write the extracted data in a csv file
        c_name = c_name.replace(" ","_")
        df.to_csv(datapath+c_name+'.csv')
        print(". done\n",end="")

    print("\n\nTime taken: {} seconds".format(round(time.time()-start, 2)))


##############################################################
## Now, we define the different functions used in the main  ##
##############################################################

# To download the tota list of countries present on worldometers
def get_country_list():

    html_page = load_data(baseurl)

    res = html_parser(html_page.content, 'a', class_='mt_a')

    country_lst = {}

    for r in tqdm(res):
        country_lst.update({r.get_text():r.get('href')})

    df_country = pd.DataFrame()

    for k,v in country_lst.items():
        df_country = df_country.append([[k,v]], ignore_index='True')
    df_country.columns = ['country','link']

    return df_country

# To load the html worldometers website
def load_data(link):
    
    try:
        html_page = requests.get(link)
    except requests.exceptions.RequestException as e:
        print (e)

    return html_page

# To parse the htmp page according to a specific tag
def html_parser(content, tag, class_=False, id_=False):
    
    bs = BeautifulSoup(content, 'html.parser')
    
    if class_:
        search = bs.find_all(tag, class_=class_)
    elif id_:
        search = bs.find_all(id=id_)
    else:
        search = bs.select(tag)
    
    return search

# To dowload the data for the selected country
def clean_data(res, info_label):

    # find the right script tag
    data = False
    for i, r in enumerate(res):
        r = str(r)
        if r.find(info_label) != -1:
            data = r.split("\n")
            break

    if data == False:
        return False, False

    # find the right line inside script tag
    for i,d in enumerate(data):
        if d.find('xAxis') != -1:
            date_range = data[i+1]
            date_range = data[i+1]
            date_range = date_range.replace("\",\"","dummy")
            date_range = date_range.replace(",","$")
            date_range = date_range.replace("dummy","\",\"")
            date_range = re.search(r"(?<=\[).*?(?=\])", date_range).group(0).split(",")
            date_range = [t.strip('\"') for t in date_range]
        elif d.find('series') != -1:
            data_range = data[i+4]
            data_range = re.search(r"(?<=\[).*?(?=\])", data_range).group(0).split(",")
            data_range = [int(t) if t != 'null' else 0 for t in data_range]

    return date_range, data_range

# execute main
if __name__ == "__main__":
    main()
