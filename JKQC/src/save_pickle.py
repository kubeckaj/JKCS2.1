def save_pickle(tosave,output_pkl,Qsplit,Qout):
  tosave = tosave.reset_index(drop=True)
  if Qsplit == 1:
    try:
      tosave.to_pickle(output_pkl)
    except:
      print("Pickle was not written down due to an error.")
    if Qout >= 1:
      if len(tosave) == 0:
        print(tosave)
        print("No files in the input!")
      else:
        print("Example output:")
        print(tosave.iloc[0])
        print("Number of files in "+output_pkl+": "+str(len(tosave)))
  else:
    if len(tosave) == 0:
      if Qout >= 1:
        print(tosave)
        print("No files in the input!")
    else:
      lengths = -(-len(tosave)//Qsplit)
      for split in range(Qsplit):
        output_pkl_split = output_pkl[:-4]+"_s"+str(split+1)+".pkl"
        start=split*lengths
        end=(split+1)*lengths
        if end > len(tosave):
          end = len(tosave)
        tosave.loc[start:end].to_pickle(output_pkl_split)
        if Qout >= 1:
          print("Number of files in "+output_pkl_split+": "+str(len(tosave.loc[start:end])))
