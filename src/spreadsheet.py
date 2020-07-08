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

# use credentials to create a client to interact with the Google Drive API
scope = ['https://spreadsheets.google.com/feeds']
creds = ServiceAccountCredentials.from_json_keyfile_name('config/client_secret.json', scope)
client = gspread.authorize(creds)

# Find a workbook by name and open the first sheet -> make sure to use the right name
# share the the spreadsheet with the credential email (beamtest2@beamtest-219608.iam.gserviceaccount.com)
sheet = client.open_by_key('1vtwJnPLbk0M1UztpSX9SZNsPYAyMCO0TnYQzD6jQWoo').sheet1


def colnum_string(n):
    string = ""
    while n > 0:
        n, remainder = divmod(n - 1, 26)
        string = chr(65 + remainder) + string
    return string


def make_desy_run_log():
    data = sheet.get_all_values()[1:]
    dic = OrderedDict()
    for row in data:
        run = row[0]
        if not run or not run.isdigit() or not row[11]:
            continue
        # row[3] += ':00' if row[3].count(':') == 1 else ''
        dic[run] = {'start': mktime(datetime.strptime('{}-{}'.format(row[1], row[2]), '%m/%d/%Y-%H:%M').timetuple()),
                    'end': mktime(datetime.strptime('{}-{}'.format(row[1], row[3]), '%m/%d/%Y-%H:%M').timetuple()),
                    'events': int(float(row[5]) * 1e6),
                    'dut0': row[7],
                    'hvsupply0': row[8],
                    'hv0': int(row[9]),
                    'current0': float(row[10]) if row[10] else '?',
                    'trim0': row[11],
                    'dut1': row[12],
                    'hvsupply1': row[13],
                    'hv1': int(row[14]),
                    'current1': float(row[15]) if row[15] else '?',
                    'trim1': row[16],
                    'angle': 0 if row[17] == '-' else int(row[17]),
                    'runplan': row[19],
                    'batch': row[20],
                    'comment': row[21]
                    }
    with open('runlog.json', 'w') as f:
        dump(dic, f, indent=2)
    print('successfully extracted data from spreadsheet and saved it as "runlog.json"')


if __name__ == '__main__':
    make_desy_run_log()
