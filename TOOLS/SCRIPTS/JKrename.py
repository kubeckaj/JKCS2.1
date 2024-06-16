def print_help():
  print("It will rename file names.")
  print("  JKrename <FILES> [OPTIONS]")
  print("OPTIONS:")
  print("  -rename X Y renames e.g. 1X2sa to 1Y2sa")
  print("EXAMPLES:")
  print("  JKrename 1sa1san.xyz -rename san sa #->2sa.xyz")
  print("  JKrename *.xyz *log -rename san sa -rename sa b")
  print("")
  print("When you are sure about the renaming, pipe it to shell to execute it:")
  print("  JKrename *.xyz *log -rename san sa -rename sa b | sh")

from os import path
from sys import argv
from os import system

QrenameWHAT = []
files = []
last = ""
for i in argv[1:]:
    #HELP
    if i == "-help" or i == "--help":
      print_help()
      exit()
    #RENAME
    if i == "-rename":
      Qmodify = 1
      Qrename = 1
      last = "-rename"
      continue
    if last == "-rename":
      remember = i
      last = "-rename2"
      continue
    if last == "-rename2":
      last = ""
      QrenameWHAT.append([remember,i])
      continue
    if path.exists(i):
      files.append(i)
      continue
    #UNKNOWN ARGUMENT
    print("I am sorry but I do not understand the argument: "+i+" [EXITING]")
    exit()

def is_nameable(input_array):
  from numpy import mod

  nameable_test = True
  if mod(len(input_array),2) == 0:
    for input_array_i in input_array[0::2]:
      if not input_array_i.isnumeric():
        nameable_test = False
        break
    for input_array_i in input_array[1::2]:
      if input_array_i.isnumeric():
        nameable_test = False
        break
  else:
    nameable_test = False
  return nameable_test

def seperate_string_number(string):
    previous_character = string[0]
    groups = []
    newword = string[0]
    for x, i in enumerate(string[1:]):
        if i == "_" or previous_character == "_":
            newword += i
        elif i.isalpha() and previous_character.isalpha():
            newword += i
        elif i.isnumeric() and previous_character.isnumeric():
            newword += i
        else:
            groups.append(newword)
            newword = i
        previous_character = i
        if x == len(string) - 2:
            groups.append(newword)
            newword = ''
    return groups

def zeros(input_array):
  output_string = ""
  skip = 0
  for i in range(len(input_array)):
    if skip == 1:
      skip = 0
      continue
    if input_array[i] == "0":
      skip = 1
      continue
    output_string += input_array[i]
  return output_string

def combine_same(input_array):
  output_array = []
  for i in range(0,len(input_array)-1,2):
    test = 0
    for j in range(0,len(output_array)-1,2):
      if input_array[i+1] == output_array[j+1]:
        test = 1
        output_array[j] = str(int(output_array[j]) + int(input_array[i]))
    if test == 0:
      output_array.append(input_array[i])
      output_array.append(input_array[i+1])
  return output_array

def replace_first_occurrence(string, old_char, new_char):
  index = string.find(old_char)  # Find the index of the first occurrence
  if index != -1:  # If the character is found
   string = string[:index] + new_char + string[index+1:]  # Replace the character
  return string

def rearrange_formula(formula):
    import re
    # Define a regular expression pattern to capture elements and their counts
    pattern = r'([A-Z][a-z]?)(\d*)'

    # Use re.findall to extract element symbols and their counts from the formula
    elements = re.findall(pattern, formula)

    # Initialize lists to store rearranged parts of the formula
    counts = []
    symbols = []

    # Iterate through the extracted elements
    for element, count in elements:
        if count == '':  # If no count is specified, set count to 1
            count = '1'
        counts.append(count)
        symbols.append(element)

    # Rearrange the elements and counts to the desired format 'count1symbol1count2symbol2...'
    rearranged_formula = ''.join(counts[i] + symbols[i] for i in range(len(elements)))

    return rearranged_formula,counts,symbols


import re
from numpy import array
for file_i in files:
  l4 = []
  file_name = file_i
  last_dot_position = file_name.rfind('.')
  file_i_BASE = file_name[0:last_dot_position]
  file_i_EXT = file_name[last_dot_position:]
  cluster_type_array = seperate_string_number(file_i_BASE.split("-")[0].split("_")[0])
  if is_nameable(cluster_type_array):
    for QrenameONE in QrenameWHAT:
      cluster_type_array = [w.replace(QrenameONE[0],QrenameONE[1]) if QrenameONE[0] == w else w for w in cluster_type_array]
    cluster_type_array = combine_same(cluster_type_array)
    l4.append("".join(cluster_type_array)+file_i_BASE.replace(file_i_BASE.split("-")[0].split("_")[0],""))
    print("mv "+file_i+" "+l4[0]+file_i_EXT)
  else:
    print("Not namable")
