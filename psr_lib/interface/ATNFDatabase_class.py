import pandas as pd
from ..utils.utils_funcs import is_float

# class ATNFDatabase:
#     def __init__(self, filename):
#         lines = []
#         with open(filename) as f:
#             lines = f.readlines()
#         dict_list = []
#         tmp_dict = {}
#         for line in lines:
#             if len(line,) == 0:
#                 print("ERROR: Zero length line")
#                 break
#             if line[0] == '@':
#                 dict_list.append(tmp_dict)
#                 tmp_dict = {}
#             else:
#                 splitted_line = line.split()
#                 if(is_float(splitted_line[1])):
#                     tmp_dict[splitted_line[0]] = float(splitted_line[1])
#                 else:
#                     tmp_dict[splitted_line[0]] = splitted_line[1]
#         self.data = pd.DataFrame(dict_list)
#         # fill unfilled but easily calculatable fields
#         self.data.loc[self.data['F0'].notna(), 'P0'] = 1 / self.data['F0']
#         self.data.loc[self.data['P0'].notna(), 'F0'] = 1 / self.data['P0']
#         self.data.loc[self.data['F0'].notna() & self.data['F1'].notna(), 'P1'] = -self.data['F1'] / self.data['F0']**2 
#         self.data.loc[self.data['P0'].notna() & self.data['P1'].notna(), 'F1'] = -self.data['P1'] / self.data['P0']**2 
#         # some additional derived parameters
#         self.data['B12'] = (self.data['P0'] * self.data['P1'] * 1e15)**0.5
#         self.data['AGE'] = self.data['P0'] / self.data['P1'] / 2
#         self.data['Q'] = self.data['P0']**(5/7) / (self.data['P1']*1e15)**(2/7)
#         self.data['L1400'] = self.data['S1400'] * self.data['DIST_DM']**2
#         self.data['L'] = 7.4e27 * self.data['L1400']
#         self.data['Edot'] = 3.95 * 1e31 * (self.data['P1'] / 1e-15)/ self.data['P0']**3
#         self.data['Eff'] = self.data['L'] / self.data['Edot']
#         # set data table index to PSR J name
#         self.data = self.data.set_index('PSRJ')


class ATNFDatabase(pd.DataFrame):
    # def __init__(self):
    #     pass
    @classmethod
    def from_ATNFdb(cls, filename): 
        lines = []
        with open(filename) as f:
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
        DataBase = cls(dict_list)
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