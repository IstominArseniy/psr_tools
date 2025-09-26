import pandas as pd
from ..utils.utils_funcs import is_float
from .psr_class import RadioPulsar

# class ATNFDatabase:
#     def __init__(self, file_name):
#         lines = []
#         with open(file_name) as f:
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


class PSRDataFrame(pd.DataFrame):
    @property
    def _constructor(self):
        return PSRDataFrame

    @property
    def _constructor_sliced(self):
        return PSRSeries
    
    @classmethod
    def from_ATNFdb(cls, file_name): 
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
    
    # @classmethod
    # def from_csv(cls, file_name):
    #     pass

    def get_psr_as_class(self, Jname, skip_warning=False):
        PSR = RadioPulsar(Jname, self.loc[Jname]['P0'], self.loc[Jname]['B12'], 45)
        if not skip_warning:
            print('No inclination angle, R, M, Ir information. Default values were assigned.')
        return PSR
    
    # def append_from_csv(self, table_name, JName_column='PSRJ'):
    #     tmp_data_frame = pd.read_csv(table_name)
    #     tmp_data_frame.set_index(JName_column, inplace=True)
    #     self.__dict__ = pd.concat([self, tmp_data_frame]).__dict__
    #     # self = pd.concat([self, tmp_data_frame]) # REDO - Do not work like that!


    

class PSRSeries(pd.Series):
    @property
    def _constructor(self):
        return PSRSeries

    @property
    def _constructor_expanddim(self):
        return PSRDataFrame
