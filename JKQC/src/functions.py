import pandas as pd

### Append to dataframe 
def df_add_append(dataframe,label,name,indexs,variables):
  if (label,name) in dataframe.columns:
    newdataframe = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
    dataframe = dataframe.append(newdataframe)
  elif dataframe.shape == (0, 0):
    dataframe = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
  else:
    dataframe[label,name] = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )

  return dataframe

### Add to dataframe via iteration
def df_add_iter(dataframe,label,name,indexs,variables):
  if (label,name) in dataframe.columns:
    df_l = dataframe.shape[0]
    var_l = len(variables)
    for i in range(var_l):
      #dataframe[label,name].iloc[df_l-var_l+i] = variables[i]
      dataframe.at[dataframe.index.values[df_l-var_l+i],(label,name)] = variables[i]
  elif dataframe.shape == (0, 0):
    dataframe = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
  else:
    dataframe[label,name] = pd.DataFrame(index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
    for i in range(len(variables)):
      #dataframe[label,name].iloc[-1-i] = variables[-1-i]
      dataframe.at[dataframe.index.values[-1-i],(label,name)] = variables[-1-i]

  return dataframe

### Append to dataframe
def df_add_insert_to_last(dataframe,label,name,indexs,variables):
  if (label,name) in dataframe.columns:
    dataframe.at[dataframe.index.values[-1],(label,name)] = variables[-1]
    #dataframe[label,name].iloc[-1] = variables[-1]
  elif dataframe.shape == (0, 0):
    dataframe = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )
  else:
    dataframe[label,name] = pd.DataFrame(variables, index = indexs, columns = pd.MultiIndex.from_tuples([(label,name)]) )

  return dataframe
