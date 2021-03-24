# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 11:42:12 2020

@author: aplomb08
"""


# Create your views here.
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

# =============================================================================

#datapath = r"data"
datapath = "/home/jnorce/Desktop/Desktop2/Covid19_minischool/worldometers/"


baseurl = "https://www.worldometers.info/coronavirus/"

start = time.time()

# for datawise data (useful for backup)
#folder_name = strftime("%Y-%b-%d", gmtime())

#try:
#    if not os.path.exists(os.path.join(datapath, folder_name)):
#        os.makedirs(os.path.join(datapath, folder_name))
#except Exception as e:
#    print("\n **error** \n")

# =============================================================================

def load_data(link):
    
    try:
        html_page = requests.get(link)
    except requests.exceptions.RequestException as e:
        print (e)

    return html_page


def html_parser(content, tag, class_=False, id_=False):
    
    bs = BeautifulSoup(content, 'html.parser')
    
    if class_:
        search = bs.find_all(tag, class_=class_)
    elif id_:
        search = bs.find_all(id=id_)
    else:
        search = bs.select(tag)
    
    return search


def clean_data(res, info_label):
    
    # find the right sciprt tag
    data = False
    for i, r in enumerate(res):
        r = str(r)
        if r.find(info_label) != -1:
            data = r.split("\n")
            # print(i)
            break
    
    if data == False:
        return False, False
    
    # find the right line inside script tag
    for i,d in enumerate(data):
#        print("xAxis: ")
#        print(d.find('xAxis'))
#        print("series: ")
#        print(d.find('series'))
#        print(data)
        if d.find('xAxis') != -1:
            date_range = data[i+1]
            date_range = date_range.replace("\",\"","dummy")
            date_range = date_range.replace(",","$")
            date_range = date_range.replace("dummy","\",\"")
#            print(date_range)
            date_range = re.search(r"(?<=\[).*?(?=\])", date_range).group(0).split("\",\"")
            date_range = [t.strip('\"') for t in date_range]
            # print(i)
        elif d.find('series') != -1:
            data_range = data[i+4]
            data_range = re.search(r"(?<=\[).*?(?=\])", data_range).group(0).split(",")
            data_range = [int(t) if t != 'null' else 0 for t in data_range]
            # print(i)
    
    return date_range, data_range

# =============================================================================
# get list of countries and links
# 1st step
# run get_country_list() and then comment it.
# then run 2nd step
# =============================================================================

# =============================================================================
# def get_country_list():
#     
#     link = "https://www.worldometers.info/coronavirus/"
#     
#     html_page = load_data(link)
#     
#     res = html_parser(html_page.content, 'a', class_='mt_a')
#     
#     country_lst = {}
#     
#     for r in tqdm(res):
#         country_lst.update({r.get_text():r.get('href')}) 
#     
#     df_country = pd.DataFrame()
#     
#     for k,v in country_lst.items():
#         df_country = df_country.append([[k,v]], ignore_index='True')
#     df_country.columns = ['country','link']
#     
#     return df_country
# 
# df_country_list = get_country_list()
# df_country_list.to_csv(os.path.join(datapath, 'country_list.csv'))

df_country_list = pd.read_csv(os.path.join(r"data", 'country_list.csv'))[['country','link']]

selected_countries = "USA ; Brazil ; UK ; Mexico ; Italy ; France ; Spain ; India ; Iran ; Peru ; Russia ; Belgium ; Germany ; Canada ; Chile ; Netherlands ; Colombia ; Sweden ; Turkey ; Pakistan ; Ecuador ; South Africa ; China ; Argentina ; Belarus ; Bolivia ; Venezuela ; Uruguay ; Paraguay ; Philippines ; Iraq ; Saudi Arabia ; Bangladesh ; Indonesia ; Ukraine ; Israel ; Qatar ; Kazakhstan ; Australia ; Lebanon ; Switzerland ; Japan ; S. Korea ;  Vietnam ; New Zealand ; Mongolia"

# =============================================================================
# get data for each country
# id, date, new_case, new_death, active, recover
# =============================================================================

path = datapath+'*.csv'
filename = [os.path.basename(x)[:-4] for x in glob.glob(path)]

# country list index that have some issue while fetching/cleaning data
country_lst_issue = []

for i, row in df_country_list.iterrows():
    
    c_name = row['country']
    link = row['link']
    if(selected_countries.find(c_name) >= 0):
        print(i, c_name, end=" ")
    else:
        continue
    
#    if c_name in filename:
#        print(i, c_name, ":)")
#        continue
#    else:
    
    
    df = pd.DataFrame()
    
    link = baseurl+link
    try:
        html_page = load_data(link)
        res = html_parser(html_page.content, 'script')
    except Exception as e:
        country_lst_issue.append([i,c_name,link])
        print(" :( :( :( ")
        continue
        
    print(".",end="")
        
#    date_range, data_daily_cases = clean_data(res, 'Daily New Cases')
#    if date_range != False:
#        df['Date'] = date_range
#        df['Daily Cases'] = data_daily_cases
#
    date_range, data_total_cases = clean_data(res, 'Total Cases')
#    print(data_total_cases)
#    print(date_range)
    if date_range != False:
        df['Date'] = date_range
        df['Total Cases'] = data_total_cases

#    date_range, data_daily_deaths = clean_data(res, 'Daily Deaths')
#    if date_range != False:
#        df['Daily Deaths'] = data_daily_deaths
#
    date_range, data_total_deaths = clean_data(res, 'Total Deaths')
    if date_range != False:
        df['Total Deaths'] = data_total_deaths
#
#    date_range, data_active = clean_data(res, 'Active Cases')
#    if date_range != False:
#        df['Active Cases'] = data_active
    
#    if not(date_range and data_daily_cases and data_daily_deaths and data_active):
#        country_lst_issue.append([i,c_name,link])
#        print(" :( :( :( ",end="")
    
    c_name = c_name.replace(" ","_")
    df.to_csv(datapath+c_name+'.csv')
    print(".",end="")
    
    print()

print("\n\nTime taken: {} seconds".format(round(time.time()-start, 2)))
