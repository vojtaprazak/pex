import urllib.request, urllib.error, urllib.parse
import json
import pickle
from numpy import *
import re, datetime
from time import strptime

def get_emdb_info(num):
    try:
        info = urllib.request.urlopen("https://www.ebi.ac.uk/pdbe/api/emdb/entry/all/EMD-"+str(num)).read()
        js_info = json.loads(info)
    except:
        js_info = 0
    return js_info

def get_all_info():
    info = {}
    for x in range(1000,10000):
        new_info = get_emdb_info(x)
        if new_info != 0:
            info.update(new_info)
    return info


class emdb_stats:

    def __init__(self, datafile):
        self.data = pickle.load(file(datafile))
        self.entry_names = array(sorted(self.data.keys()))
        self.res_stat = "processing:reconstruction:resolution_by_author"
        self.date_stat = "deposition:deposition_date"


    def get_all_entry_names(self):
        return self.entry_names


    def get_stats(self, stat, key_list=False):
        """stat is the relevant keys for the stat wanted, separated by colons (:)"""
        stat = stat.split(':')
        stat_list = []
        fail_list = []
        fail = False
        entry_names = self.get_all_entry_names()
        if key_list != False:
            entry_names = key_list
        for ent in entry_names:
            new_stat = self.data[ent][0]
            for x in range(len(stat)):
                try:
                    if stat[x] == '0':
                        new_stat = new_stat[0]
                    else:
                        new_stat = new_stat[stat[x]]
                except KeyError:
                    fail_list.append(ent)
                    fail = True
            if not fail:
                stat_list.append(array([ent, new_stat]))
            fail = False
        return array(stat_list)

    def get_entries_by_stat_value(self, stat, stat_value, key_list=False):
        all_stat_vals = self.get_stats(stat, key_list=key_list)
        #picked_stat_vals = array([x for x in all_stat_vals if type(stat_value)(x[1]) == stat_value])
        picked_stat_vals = array([x for x in all_stat_vals if re.search(stat_value, x[1])])
        return picked_stat_vals


    def get_unique_values_of_stat(self, stat, key_list=False):
        all_stats = self.get_stats(stat, key_list=key_list)
        if len(all_stats) != 0:
            return unique(all_stats[:,1])
        else:
            return array([])


    def get_two_stats(self, stat1, stat2, key_list=False):
        stat1_list = self.get_stats(stat1, key_list=key_list)
        stat2_list = self.get_stats(stat2, key_list=stat1_list[:,0])
        stat1_list = self.get_stats(stat1, key_list=stat2_list[:,0])
        return vstack((stat1_list[:,0], stat1_list[:,1], stat2_list[:,1]))


    def compare_to_date(self, stat1, key_list=False):
        a = self.get_two_stats(self.date_stat, stat1, key_list=key_list)
        b = [strptime(x, '%Y-%m-%d') for x in a[1]]
        return a[0], b, a[2]

    def compare_to_res(self, stat1, key_list=False):
        a = self.get_two_stats(self.res_stat, stat1, key_list=key_list)
        b = [float(x) for x in a[1]]
        return a[0], b, a[2]


        
