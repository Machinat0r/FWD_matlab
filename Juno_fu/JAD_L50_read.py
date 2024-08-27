#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 00:52:44 2023

------- written by Wending Fu in Beijing ------
"""
import re
import sys
import scipy
import struct

def extract_lbl(lbl_file_path):
    with open(lbl_file_path, "r") as lbl_file:
        lbl_content = lbl_file.read()
        
        # extract filename
        dat_filename_match = re.search(r'^\s*\^TABLE\s*=\s*"([^"]+\.DAT)"', lbl_content, re.MULTILINE)
        if dat_filename_match:
            dat_filename = dat_filename_match.group(1)
        
        # extract record info
        bytes_per_record_match = re.search(r'RJW,\s*BYTES_PER_RECORD,\s*(\d+)', lbl_content)
        objects_per_record_match = re.search(r'RJW,\s*OBJECTS_PER_RECORD,\s*(\d+)', lbl_content)
        if bytes_per_record_match and objects_per_record_match:
            bytes_per_record = int(bytes_per_record_match.group(1))
            objects_per_record = int(objects_per_record_match.group(1))
        
        # extract object descriptions
        object_pattern = re.compile(r'(OBJECT\s*=\s*CONTAINER.+?|OBJECT\s*=\s*COLUMN.+?)END_OBJECT', re.MULTILINE | re.DOTALL)
        object_matches = object_pattern.findall(lbl_content)
        object_descriptions = {}
        for object_match in object_matches:
            object_name_match = re.search(r'NAME\s*=\s*([^\s]+)', object_match)
            if object_name_match:
                object_name = object_name_match.group(1)
                object_descriptions[object_name] = object_match.strip()
                    
        return dat_filename, bytes_per_record, objects_per_record, object_descriptions


def read_dat_file(dat_file_path, bytes_per_record):
    with open(dat_file_path, "rb") as dat_file:
        records = []
        while True:
            record_data = dat_file.read(bytes_per_record)
            if not record_data:
                break
            records.append(record_data)
    return records

def unpack_object(obj_description, data):
    object_name = re.findall(r'NAME\s*=\s*([^\s]+)', obj_description)[-1]
    start_byte = int(re.search(r'START_BYTE\s*=\s*(\d+)', obj_description).group(1))
    if object_name == 'ISSUES_BITS':
        bytes_length = int(re.findall(r'BYTES\s*=\s*(\d+)', obj_description)[0])
        data_type = 'I'
        data_length2 = 1
        data_length1 = '1'
    else:
        bytes_length = int(re.findall(r'BYTES\s*=\s*(\d+)', obj_description)[-1])
        data_type = re.search('RJW,\s*' + object_name + ',\s*([A-Za-z]+),\s', obj_description).group(1)
        data_length1 = re.search('RJW,\s*' + object_name + ',\s*.+\s([0-9]+)\s\*/', obj_description).group(1)
        data_length2 = int(re.search('RJW,\s*' + object_name + ',\s*.+\s([0-9]+),\s*[0-9]+\s\*/', obj_description).group(1))
        
    format_string = f"{data_length1}{data_type}"
    if data_length2 == 1:
        unpacked_data = struct.unpack(format_string, data[start_byte - 1:start_byte + bytes_length - 1])
    else:
        unpacked_data = []
        for row in range(data_length2):
            unpacked_data.append(struct.unpack(format_string, data[start_byte + row * int(data_length1) * bytes_length - 1:start_byte + (row+1) * int(data_length1) * bytes_length - 1]))
            
    return unpacked_data

def unpack_records(records, object_descriptions, var_name):
    unpacked_records = []

    for record in records:
        if len(record) != bytes_per_record:
            Warning('Wrong bytes per record. Please check the Dat file.')
                
        record_data = {}
        if '*' not in var_name:
            for object_name in var_name:
                try:
                    obj_description = object_descriptions[object_name]
                    unpacked_data = unpack_object(obj_description, record)
                    record_data[object_name] = unpacked_data
                except:
                    Warning('Wrong Var Name!!!')
        else:
            for object_name, obj_description in object_descriptions.items():
                unpacked_data = unpack_object(obj_description, record)
                record_data[object_name] = unpacked_data
        unpacked_records.append(record_data)

    return unpacked_records

def organize_data(unpacked_records):
    organized_data = {key: [] for key in unpacked_records[0]}
    for item in unpacked_records:
        for key, value in item.items():
            organized_data[key].append(value)
    
    return organized_data


if __name__ == "__main__":
    #lbl_file_path = "/Users/fwd/Documents/MATLAB/Code/wcq/20170319_Rolling-pin/JADE/JAD_L50_LRS_ELC_ANY_DEF_2019307_V01.LBL"
    #var_name = 'DIM0_UTC&MAG_VECTOR_DESPUN'.split('&')
    lbl_file_path = sys.argv[1]
    var_name = sys.argv[2].split('&')
    temp_dir = '/'.join(lbl_file_path.split('/')[:-1]) + '/'
    
    dat_filename, bytes_per_record, objects_per_record, object_descriptions= extract_lbl(lbl_file_path)
    if objects_per_record != len(object_descriptions):
        Warning('Wrong objects number. Please check the LBL file.')
        
    dat_file_path = temp_dir + dat_filename  # .DAT & .LBL should be in the same dir
    records = read_dat_file(dat_file_path, bytes_per_record)

    unpacked_records = unpack_records(records, object_descriptions, var_name)
    organized_data = organize_data(unpacked_records)
    
    
    scipy.io.savemat(temp_dir + 'organized_data.mat', organized_data);
