import pandas as pd
from ..utils.utils_funcs import is_float
from .psr_class import RadioPulsar

def read_from_ATNF_DB(file_name): 
        lines = []
        with open(file_name) as f:
            lines = f.readlines()
        dict_list = []
        tmp_dict = {}
        for line in lines:
            if len(line,) == 0:
                print("ERROR: Zero length line")
                break
            if line[0] == '@':
                dict_list.append(tmp_dict)
                tmp_dict = {}
            else:
                splitted_line = line.split()
                if(is_float(splitted_line[1])):
                    tmp_dict[splitted_line[0]] = float(splitted_line[1])
                else:
                    tmp_dict[splitted_line[0]] = splitted_line[1]
        DataBase = pd.DataFrame(dict_list)
        # fill unfilled but easily calculatable fields
        DataBase.loc[DataBase['F0'].notna(), 'P0'] = 1 / DataBase['F0']
        DataBase.loc[DataBase['P0'].notna(), 'F0'] = 1 / DataBase['P0']
        DataBase.loc[DataBase['F0'].notna() & DataBase['F1'].notna(), 'P1'] = -DataBase['F1'] / DataBase['F0']**2 
        DataBase.loc[DataBase['P0'].notna() & DataBase['P1'].notna(), 'F1'] = -DataBase['P1'] / DataBase['P0']**2 
        # some additional derived parameters
        DataBase['B12'] = (DataBase['P0'] * DataBase['P1'] * 1e15)**0.5
        DataBase['AGE'] = DataBase['P0'] / DataBase['P1'] / 2
        DataBase['Q'] = DataBase['P0']**(5/7) / (DataBase['P1']*1e15)**(2/7)
        DataBase['L1400'] = DataBase['S1400'] * DataBase['DIST_DM']**2
        DataBase['L'] = 7.4e27 * DataBase['L1400']
        DataBase['Edot'] = 3.95 * 1e31 * (DataBase['P1'] / 1e-15)/ DataBase['P0']**3
        DataBase['Eff'] = DataBase['L'] / DataBase['Edot']
        # set data table index to PSR J name
        DataBase = DataBase.set_index('PSRJ')
        return DataBase

def read_from_csv(file_name, JName_column='PSRJ'):
    DataFrame = pd.read_csv(file_name)
    DataFrame.set_index(JName_column, inplace=True)
    return DataFrame

def get_psr_from_table(self, Jname, skip_warning=False):
    PSR = RadioPulsar(Jname, self.loc[Jname]['P0'], self.loc[Jname]['B12'], 45)
    if not skip_warning:
        print('No inclination angle, R, M, Ir information. Default values were assigned.')
    return PSR