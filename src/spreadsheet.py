#!/usr/bin/env python
# --------------------------------------------------------
#       small script to read data from google spreadsheet
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

import gspread
from oauth2client.service_account import ServiceAccountCredentials
from collections import OrderedDict
from datetime import datetime
from json import dump
from time import mktime

year = '2020'

# use credentials to create a client to interact with the Google Drive API
scope = ['https://spreadsheets.google.com/feeds']
creds = ServiceAccountCredentials.from_json_keyfile_name('config/client_secret.json', scope)
client = gspread.authorize(creds)

# Find a workbook by name and open the first sheet -> make sure to use the right name
# share the the spreadsheet with the credential email (beamtest2@beamtest-219608.iam.gserviceaccount.com)
sheet = client.open_by_key('1vtwJnPLbk0M1UztpSX9SZNsPYAyMCO0TnYQzD6jQWoo').sheet1

# Extract and print all of the values
# list_of_hashes = sheet.get_all_records()
# print(list_of_hashes)


def colnum_string(n):
    string = ""
    while n > 0:
        n, remainder = divmod(n - 1, 26)
        string = chr(65 + remainder) + string
    return string


col2num = lambda col: reduce(lambda x, y: x * 26 + y, [ord(c.upper()) - ord('A') + 1 for c in col])


def make_desy_run_log():
    data = sheet.get_all_values()[1:]
    dic = OrderedDict()
    for row in data:
        run = row[0]
        if not run or not run.isdigit() or not row[11]:
            continue
        row[3] += ':00' if row[3].count(':') == 1 else ''
        dic[run] = {'start': mktime(datetime.strptime('{}{}{}'.format(year, row[1], row[2]), '%Y%a, %b %d %H:%M').timetuple()),
                    'stop': mktime(datetime.strptime('{}{}{}'.format(year, row[1], row[3]), '%Y%a, %b %d %H:%M:%S').timetuple()),
                    'events': int(float(row[5]) * 1e6),
                    'dut1': row[7],
                    'hvname1': row[8],
                    'hv1': int(row[9]),
                    'current1': float(row[10]) if row[10] else '?',
                    'dut2': row[11],
                    'hvname2': row[12],
                    'hv2': int(row[13]),
                    'current2': float(row[14]) if row[14] else '?',
                    'angle': 0 if row[15] == '-' or '/' in row[15] else int(row[15]),
                    'runplan': row[17],
                    'batch': row[18],
                    'comment': row[19]
                    }
    with open('runlog.json', 'w') as f:
        dump(dic, f, indent=2)
