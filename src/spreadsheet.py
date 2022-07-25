#!/usr/bin/env python
# --------------------------------------------------------
#       small script to read data from Google spreadsheet
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

import gspread
from oauth2client.service_account import ServiceAccountCredentials
from datetime import datetime
from json import dump
from time import mktime
from utility.utils import Dir, array

# use credentials to create a client to interact with the Google Drive API
scope = ['https://spreadsheets.google.com/feeds']
creds = ServiceAccountCredentials.from_json_keyfile_name(Dir.joinpath('config', 'client_secret.json'), scope)
client = gspread.authorize(creds)

# Find a workbook by name and open the first sheet -> make sure to use the right name
# share the spreadsheet with the credential email (beamtest2@beamtest-219608.iam.gserviceaccount.com)


def colnum_string(n):
    string = ""
    while n > 0:
        n, remainder = divmod(n - 1, 26)
        string = chr(65 + remainder) + string
    return string


def make_timestamp(date, time):
    return int(mktime(datetime.strptime(f'{date}{time}', '%m/%d/%Y%H:%M').timetuple()))


def make_desy_run_log():
    sheet = client.open_by_key('1vtwJnPLbk0M1UztpSX9SZNsPYAyMCO0TnYQzD6jQWoo').sheet1
    data = sheet.get_all_values()[1:]
    dic = {}
    for row in data:
        run = row[0]
        if not run or not run.isdigit() or not row[11]:
            continue
        # row[3] += ':00' if row[3].count(':') == 1 else ''
        dut_data = array(row[7:17]).reshape((2, -1)).T.tolist()  # two DUTs
        dut_dict = {name: d for name, d in zip(['duts', 'hv supplies', 'hv', 'current', 'trim'], dut_data)}
        dic[run] = {'start': make_timestamp(row[1], row[2]),
                    'end': make_timestamp(row[1], row[3]),
                    'events': int(float(row[5]) * 1e6),
                    **dut_dict,
                    'angle': 0 if row[17] == '-' else int(row[17]),
                    'runplan': row[19],
                    'batch': row[20],
                    'comment': row[21]
                    }
    with open('runlog.json', 'w') as f:
        dump(dic, f, indent=2)
    print('successfully extracted data from spreadsheet and saved it as "runlog.json"')


info = {'2018-09': {'key': '1KoDi9OLU0SiqtLvTGgglQGHm51xn5J0ErAs89h6a-Rc',
                    'hv': {'II6-A2': '1-5', 'CMS04': '1-4', 'Si-D7': '1-4'},
                    'n': [3, 4, 5, 28]},  # [row0 for dut, n dut, n dut rows, row0 for info]
        '2018-10': {'key': '1t-MXNW0eN9tkGZSakfPdmnd_wcq4cX14Nw0bQ2ma_OQ',
                    'hv': {'II6-A2': '2-1', 'CMS04': '2-3', 'Si-D8': '2-2', 'II6-B6': '2-3'},
                    'n': [4, 4, 5, 29]}
        }


def make_cern_run_log(tc='2018-10'):
    sheet = client.open_by_key(info[tc]['key']).worksheet('KARTEL')
    data = sheet.get_all_values()[1:]
    hv = info[tc]['hv']
    exclude = ['FEI4', '1x5']
    dic = {}
    for row in data:
        run = row[1]
        r0, ndut, nrow, r1 = info[tc]['n']
        if not run or not run.isdigit() or not row[r1 + 1] or not row[r1 + 2] or not row[0].isdigit():
            continue
        dut_data = array(row[r0:r0 + ndut * nrow]).reshape((ndut, -1))  # up to four DUTs
        dut_data = array([[w.strip(' ') for w in lst] for lst in dut_data if lst[0] and not any([w in lst[0] for w in exclude])]).T.tolist()  # remove empty and excluded DUTs
        dut_dict = {name: d for name, d in zip(['duts', 'hv', 'current', 'temp', 'angle'], dut_data)}
        if not dut_dict:
            continue
        dic[run] = {'telescope run': int(row[0]),
                    'start': make_timestamp(row[r1 + 1], row[r1 + 2]),
                    'end': make_timestamp(row[r1 + 1], row[r1 + 3]),
                    'events': int(row[r1 + 6]) if row[r1 + 6] else 0,
                    **dut_dict,
                    'hv supplies': [hv[dut] if dut in hv else '' for dut in dut_dict['duts']],
                    'status': row[r1 + 5],
                    'batch': row[r1],
                    'comment': row[r1 + 8],
                    }
    with open('runlog.json', 'w') as f:
        dump(dic, f, indent=2)
    print('successfully extracted data from spreadsheet and saved it as "runlog.json"')


if __name__ == '__main__':
    pass
    # make_desy_run_log()
